#include "kmerops.hpp"
#include "logger.hpp"
#include "dnaseq.hpp"
#include "timer.hpp"
#include "paradissort.hpp"
#include <cstring>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <upcxx/upcxx.hpp>
#include <thread>

std::vector<upcxx::persona>* g_personas;
KmerSeedBuckets* g_bucket;
thread_local int done_cnt = 0;
std::atomic<bool> g_send_done[THREADCNT_SEND];
std::atomic<bool> g_recv_done[THREADCNT_RECV];


void store_into_local(const int& tid, const size_t& buf_sz, const upcxx::view<uint8_t>& buf_in_rpc, const bool last_send){
    const uint8_t* addr = buf_in_rpc.data();
    const uint8_t* end = addr + buf_sz;

    // std::cout<<"Storing into "<< tid << " with size " << buf_sz << " last_send? "<<last_send << std::endl;
    
    while(addr < end){
        TKmer kmer((void*)addr);                
        ReadId readid = *((ReadId*)(addr + TKmer::NBYTES));
        PosInRead pos = *((PosInRead*)(addr + TKmer::NBYTES + sizeof(ReadId)));
        (*g_bucket)[tid].emplace_back(kmer, readid, pos);
        addr += TKmer::NBYTES + sizeof(ReadId) + sizeof(PosInRead);
    }

    if (last_send){
        done_cnt++;
    }
}

void KmerParserHandler::encode(std::vector<uint8_t>& t, const TKmer& kmer, PosInRead posinread, ReadId readid, size_t& buf_sz) {

    uint8_t* addr = &t[buf_sz];
    memcpy(addr, kmer.GetBytes(), TKmer::NBYTES);
    memcpy(addr + TKmer::NBYTES, &readid, sizeof(ReadId));
    memcpy(addr + TKmer::NBYTES + sizeof(ReadId), &posinread, sizeof(PosInRead));
    buf_sz += (sizeof(TKmer) + sizeof(ReadId) + sizeof(PosInRead));

}

upcxx::future<> KmerParserHandler::send(int dst, bool last_send) {
        
        int dst_process = dst / ntasks;
        // std::cout<<"From process:"<<upcxx::rank_me()<<" sending to process:"<<dst_process<< " local taskid:"<< dst % ntasks << " last_send:" << last_send << std::endl;

        auto task_done = upcxx::rpc(
            dst_process,
            [](const int& taskid, const size_t& buf_sz, const upcxx::view<uint8_t>& buf_in_rpc, bool last_send) {
                upcxx::persona& persona = (*g_personas)[taskid % THREADCNT_RECV];

                auto task_done = persona.lpc([taskid, buf_sz, &buf_in_rpc, last_send](){
                    store_into_local(taskid, buf_sz, buf_in_rpc, last_send);
                });
                return task_done;   /* returning this future ensures that the lifetime of buffer is extended*/
            },
            dst % ntasks, buf_sz[dst], upcxx::make_view(&sendbufs[dst][0], &sendbufs[dst][buf_sz[dst]]), last_send);

        buf_sz[dst] = 0;
        return task_done;
}

void KmerParserHandler::wait() {
    for (int dst = 0; dst < nprocs * ntasks; dst++) {
        futures = upcxx::when_all(futures, send(dst, true));
    }
    futures.wait();
}


void KmerParserHandler::operator()(const TKmer& kmer, size_t kid, size_t rid) {
    int dst = GetKmerOwner(kmer, nprocs * ntasks);
    encode(sendbufs[dst], kmer, static_cast<PosInRead>(kid), static_cast<ReadId>(rid) + readoffset, buf_sz[dst]);

    if (buf_sz[dst] >= batch_sz) {
        futures = when_all(futures, send(dst));
    }
}

std::unique_ptr<KmerSeedBuckets>
exchange_kmer(const DnaBuffer& myreads,
     int thr_per_task)
{

    #if LOG_LEVEL >= 3
    Timer timer();
    #endif

    int myrank = upcxx::rank_me();
    int nprocs = upcxx::rank_n();
    bool single_node = (nprocs == 1);

    Logger logger;
    std::ostringstream rootlog;

    /* parallel settings */
    
    int ntasks = omp_get_max_threads() / thr_per_task;
    if (ntasks < 1) {
        ntasks = 1;
        thr_per_task = omp_get_max_threads();
    }
    
    int nthr_send = THREADCNT_SEND;
    int nthr_recv = std::min(THREADCNT_RECV, ntasks);

    if (ntasks % nthr_recv != 0) {
        if (myrank==0){
            logger() << "Warning: ntasks is not divisible by nthr_recv. work may be imbalanced.";
            logger.flush(0);
        }
    }


    logger() << "ntasks: " << ntasks << " nthr_send: " << nthr_send << " nthr_recv: " << nthr_recv;
    logger.flush("Thread Info");


    size_t numreads = myreads.size();     /* Number of locally stored reads */

    #if LOG_LEVEL >= 3
    timer.start();
    #endif

    size_t readoffset = 0;
    
    upcxx::dist_object<size_t> readoffset_dist(numreads);
    for(int i = 0; i < myrank; i++)
        readoffset += readoffset_dist.fetch(i).wait();
    upcxx::barrier();

    std::vector<upcxx::persona> personas(nthr_recv);

    std::vector<std::thread> send_threads;
    std::vector<std::thread> recv_threads;

    g_personas = &personas;
    KmerSeedBuckets* bucket = new KmerSeedBuckets(ntasks);
    // maybe reserve space for the bucket will help
    g_bucket = bucket;

    upcxx::barrier();

    /* launch the receiving threads first */
    for (int i = 0; i < nthr_recv; i++ ){
        g_recv_done[i] = false;
        int my_taskcnt = ntasks / nthr_recv;
        if (i < ntasks % nthr_recv) my_taskcnt++;
        int expected_done = my_taskcnt * nprocs * nthr_send; 

        logger() <<"Launching recv thread "<< i << " with expected done "<< expected_done << std::endl;
        logger.flush("Launching recv threads");

        recv_threads.emplace_back(std::thread([](int expected_done, upcxx::persona& persona, int i){

            upcxx::persona_scope scope(persona);
            while(expected_done != done_cnt)    /* done_cnt is a thread local var */
                upcxx::progress();

            g_recv_done[i] = true;
            upcxx::discharge();

        }, expected_done, std::ref(personas[i]), i));
    }


    for (int i = 0; i < nthr_send; i++ ){
        g_send_done[i] = false;
        size_t st = numreads / nthr_send * i;
        size_t ed = numreads / nthr_send * (i + 1);
        if (i == nthr_send - 1) ed = numreads;

        send_threads.emplace_back(std::thread([&](size_t st, size_t ed, int i, const DnaBuffer& myreads){

            KmerParserHandler handler(ntasks, nprocs, BATCH_SIZE, readoffset);
            ForeachKmer(myreads, handler, st, ed);
            handler.wait();

            // std::cout<<"Finished! Send thr no:"<<i<< " process rank: "<<myrank<< std::endl;
            g_send_done[i] = true;
            upcxx::discharge();

        }, st, ed, i, std::ref(myreads)));

    }

    while (true) {
        upcxx::progress();

        bool all_done = true;
        for (int i = 0; i < nthr_send; i++) {
            if (g_send_done[i] == false) {
                all_done = false;
                break;
            }
        }
        for(int i = 0; i < nthr_recv; i++) {
            if (g_recv_done[i] == false) {
                all_done = false;
                break;
            }
        }

        if (all_done == true) {
            for (int i = 0; i < nthr_send; i++) {
                send_threads[i].join();
            }
            for (int i = 0; i < nthr_recv; i++) {
                recv_threads[i].join();
            }
            break;
        }
    }

    return std::unique_ptr<KmerSeedBuckets>(bucket);
}


std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, int thr_per_task)
{
    int myrank = upcxx::rank_me();
    int nprocs = upcxx::rank_n();

    Logger logger;
    std::ostringstream rootlog;

    int ntasks = omp_get_max_threads() / thr_per_task;
    if (ntasks < 1){
        ntasks = 1;
        thr_per_task = omp_get_max_threads();
    }

    assert(ntasks == recv_kmerseeds->size());

    omp_set_nested(1);
    omp_set_num_threads(ntasks);

    uint64_t task_seedcnt[ntasks];
    uint64_t valid_kmer[ntasks];
    KmerList kmerlists[ntasks];

    #if LOG_LEVEL >= 2
    logger() <<"task num: "<< ntasks << " \tthread per task: " << thr_per_task <<" \ttask size: ";
    for (int i = 0; i < ntasks; i++) {
        task_seedcnt[i] = (*recv_kmerseeds)[i].size();
        logger() << task_seedcnt[i] << " ";
    }
    logger.flush("Parallel tasks for kmer sorting and counting:");
    #endif

    #if LOG_LEVEL >= 3
    timer.start();
    #endif

    #pragma omp parallel 
    {
        int tid = omp_get_thread_num();
        paradis::sort<KmerSeedStruct, TKmer::NBYTES>((*recv_kmerseeds)[tid].data(), (*recv_kmerseeds)[tid].data() + task_seedcnt[tid], thr_per_task);
    

    #if LOG_LEVEL >= 3
    }
    timer.stop_and_log("Shared memory parallel K-mer sorting");
    timer.start();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
    #endif


        kmerlists[tid].reserve(uint64_t(task_seedcnt[tid] / LOWER_KMER_FREQ));    // This should be enough
        TKmer last_mer = (*recv_kmerseeds)[tid][0].kmer;
        uint64_t cur_kmer_cnt = 1;
        valid_kmer[tid] = 0;

        for(size_t idx = 1; idx < task_seedcnt[tid]; idx++) {
            TKmer cur_mer = (*recv_kmerseeds)[tid][idx].kmer;
            if (cur_mer == last_mer) {
                cur_kmer_cnt++;
            } else {
                if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {

                    kmerlists[tid].push_back(KmerListEntry());
                    KmerListEntry& entry    = kmerlists[tid].back();

                    TKmer& kmer             = std::get<0>(entry);
                    READIDS& readids        = std::get<1>(entry);
                    POSITIONS& positions    = std::get<2>(entry);
                    int& count              = std::get<3>(entry);

                    count = cur_kmer_cnt;
                    kmer = last_mer;

                    for (size_t j = idx - count; j < idx; j++) {
                        readids[j - idx + count] = (*recv_kmerseeds)[tid][j].readid;
                        positions[j - idx + count] = (*recv_kmerseeds)[tid][j].posinread;
                    }

                    valid_kmer[tid]++;
                }

                cur_kmer_cnt = 1;
                last_mer = cur_mer;
            }
        }

        // deal with the last kmer 
        if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {
            kmerlists[tid].push_back(KmerListEntry());
            KmerListEntry& entry         = kmerlists[tid].back();

            TKmer& kmer             = std::get<0>(entry);
            READIDS& readids        = std::get<1>(entry);
            POSITIONS& positions    = std::get<2>(entry);
            int& count              = std::get<3>(entry);

            count = cur_kmer_cnt;
            kmer = last_mer;

            for (size_t j = task_seedcnt[tid] - count; j < task_seedcnt[tid]; j++) {
                readids[j - task_seedcnt[tid] + count] = (*recv_kmerseeds)[tid][j].readid;
                positions[j - task_seedcnt[tid] + count] = (*recv_kmerseeds)[tid][j].posinread;
            }

            valid_kmer[tid]++;
        }
    }

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer counting");
    timer.start();
    #endif

    #if LOG_LEVEL >= 2
    for (int i = 0; i < ntasks; i++) {
        logger() << valid_kmer[i] << " ";
    }
    logger.flush("Valid kmer from tasks:");
    #endif


    uint64_t valid_kmer_total = std::accumulate(valid_kmer, valid_kmer + ntasks, (uint64_t)0);
    KmerList* kmerlist = new KmerList();
    kmerlist->reserve(valid_kmer_total);

    // yfli: actually we don't need to copy it into one single kmerlist
    for (int i = 0; i < ntasks; i++) {
        kmerlist->insert(kmerlist->end(), kmerlists[i].begin(), kmerlists[i].end());
    }

    logger() << valid_kmer_total;
    logger.flush("Valid kmer for process:");

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer copying");
    #endif

    return std::unique_ptr<KmerList>(kmerlist);
}

int GetKmerOwner(const TKmer& kmer, int nprocs) {
    uint64_t myhash = kmer.GetHash();
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / std::numeric_limits<uint64_t>::max();
    assert(owner >= 0 && owner < static_cast<int>(nprocs));
    return static_cast<int>(owner);
}

void print_kmer_histogram(const KmerList& kmerlist) {
    #if LOG_LEVEL >= 2

    Logger logger;
    int maxcount = std::accumulate(kmerlist.cbegin(), kmerlist.cend(), 0, [](int cur, const auto& entry) { return std::max(cur, std::get<3>(entry)); });

    maxcount = upcxx::reduce_all(maxcount, upcxx::op_fast_max).wait();

    std::vector<int> histo(maxcount+1, 0);
    std::vector<int> g_histo(maxcount+1, 0);

    for(size_t i = 0; i < kmerlist.size(); ++i)
    {
        int cnt = std::get<3>(kmerlist[i]);
        assert(cnt >= 1);
        histo[cnt]++;
    }

    upcxx::reduce_all(histo.data(), g_histo.data(), maxcount+1, upcxx::op_fast_add).wait();

    int myrank = upcxx::rank_me();

    if (!myrank)
    {
        std::cout << "#count\tnumkmers" << std::endl;

        for (int i = 1; i < g_histo.size(); ++i)
        {
            if (g_histo[i] > 0)
            {
                std::cout << i << "\t" << g_histo[i] << std::endl;
            }
        }
        std::cout << std::endl;
    }

    upcxx::barrier();
    #endif
}
