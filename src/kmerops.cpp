#include "kmerops.hpp"
#include "logger.hpp"
#include "dnaseq.hpp"
#include "timer.hpp"
#include "paradissort.hpp"
#include "memcheck.hpp"
#include "supermer.hpp"
#include "compiletime.h"
#include <cstring>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <thread>
#include <deque>


bool BatchSender::write_sendbuf(uint8_t* addr, int procid) {

    size_t task = loc_task[procid];
    size_t thr = loc_thr[procid];
    size_t idx = loc_idx[procid];

    if (task >= ntasks) return true;    /* all data have already been sent. */

    std::vector<KmerSeedStruct> *v = &(kmerseeds[thr][procid * ntasks + task]);
    size_t max_idx = v->size();
    uint8_t* addr_limit = addr + send_limit;    /* No more kmers should be written if going past this location */

    while (addr <= addr_limit && task < ntasks) {

        encode_kmerseed_basic(addr, (*v)[idx]);
        idx++;

        /* Carry over to the next group to be sent */
        if (idx >= max_idx) {
            idx = 0;
            thr++;
            if (thr >= nthrs) {
                thr = 0;
                task++;
            }
            if (task >= ntasks) break;

            v = &(kmerseeds[thr][procid * ntasks + task]);
            max_idx = v->size();
        }
    }

    /* update the location of the next sending item */
    loc_idx[procid] = idx;
    loc_thr[procid] = thr;
    loc_task[procid] = task;

    return (task >= ntasks);
}


void BatchSender::write_sendbufs(uint8_t* addr) {
    bool finished = true;

    #pragma omp parallel for num_threads(2)
    for(int i = 0; i < nprocs; i++) {
        bool myfinished = write_sendbuf(addr + batch_sz * i, i);
        
        #pragma omp critical
        {
            if (!myfinished) finished = false;
        }
    }

    /* write the last byte of each buf, which indicates the sending state of the current process */
    if (finished) {
        for (int i = 0; i < nprocs; i++)
            addr[batch_sz * i + batch_sz - 1] = 1;
    } else {
        for (int i = 0; i < nprocs; i++)
            addr[batch_sz * i + batch_sz - 1] = 0;
    }
}


uint8_t* BatchSender::progress() {
    if (status != BATCH_SENDING) return 0;

    #if LOG_LEVEL >= 3
    Timer timer(comm);
    Logger logger(comm);
    timer.start();
    double elapsed;
    #endif

    /* write the sending buffer for the current round */
    write_sendbufs(sendtmp_y);

    #if LOG_LEVEL >= 3
    elapsed = timer.get_elapsed();
    timer.stop_and_log("write_sendbufs");
    logger() <<elapsed;
    logger.flush("write_sendbufs");
    timer.start();
    #endif

    /* complete the sending and receiving of data for the current round */
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    #if LOG_LEVEL >= 3
    elapsed = timer.get_elapsed();
    timer.stop_and_log("MPI_Wait");
    logger() << elapsed;
    logger.flush("MPI_Wait");
    timer.start();
    #endif



    std::swap(sendtmp_x, sendtmp_y);
    std::swap(recvtmp_x, recvtmp_y); 

    /* check if this process has sent all its kmers */
    bool finished = true;
    for (int i = 0; i < nprocs; i++)
        if (recvtmp_y[batch_sz * i + batch_sz - 1] == false) {
            finished = false;
            break;
        }

    if (finished) {
        status = BATCH_DONE;
        return recvtmp_y;
    }

    /* start sending and receiving for the current round */
    MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);
    return recvtmp_y;
}

uint8_t* BatchSender::progress_basic_overlap() {
    if (status != BATCH_SENDING) return 0;

    #if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
    #endif

    /* complete the sending and receiving of data for the current round */
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("MPI_Wait");
    timer.start();
    #endif

    /* write the sending buffer for the current round */
    write_sendbufs(sendtmp_y);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("write_sendbufs");
    timer.start();
    #endif


    std::swap(sendtmp_x, sendtmp_y);
    std::swap(recvtmp_x, recvtmp_y); 

    /* check if this process has sent all its kmers */
    bool finished = true;
    for (int i = 0; i < nprocs; i++)
        if (recvtmp_y[batch_sz * i + batch_sz - 1] == false) {
            finished = false;
            break;
        }

    if (finished) {
        status = BATCH_DONE;
        return recvtmp_y;
    }

    /* start sending and receiving for the current round */
    MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);
    return recvtmp_y;
}

uint8_t* BatchSender::progress_basic_a() {
    if (status != BATCH_SENDING) return 0;

    #if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
    #endif

    /* complete the sending and receiving of data for the current round */
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("MPI_Wait");
    #endif

    return recvtmp_x;
}


void BatchSender::progress_basic_b() {
    #if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
    #endif

    /* write the sending buffer for the current round */
    write_sendbufs(sendtmp_y);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("write_sendbufs");
    timer.start();
    #endif


    std::swap(sendtmp_x, sendtmp_y);
    std::swap(recvtmp_x, recvtmp_y); 

    /* check if this process has sent all its kmers */
    bool finished = true;
    for (int i = 0; i < nprocs; i++)
        if (recvtmp_y[batch_sz * i + batch_sz - 1] == false) {
            finished = false;
            break;
        }

    if (finished) {
        status = BATCH_DONE;
        return;
    }

    /* start sending and receiving for the current round */
    MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);
    return;
}


void BatchReceiver::receive(uint8_t* recvbuf){
    #if LOG_LEVEL >= 3
    Timer timer(MPI_COMM_WORLD);
    timer.start();
    #endif

    for(int i = 0; i < nprocs; i++) 
    {
        if (recv_task_id[i] >= taskcnt) continue;                       /* all data for this task has been received */

        /* current receiving information */
        int working_task = recv_task_id[i];             
        size_t cnt = recvcnt[i];
        size_t max_cnt = expected_recvcnt[ taskcnt * i + working_task ];
        uint8_t* addrs2read = recvbuf + batch_size * i;
        uint8_t* addr_limit = addrs2read + send_limit;

        while(addrs2read <= addr_limit && working_task < taskcnt) {
            decode_basic(addrs2read, kmerseeds[working_task]);
            cnt++;

            /* Carry over to the next receiving group*/
            if (cnt >= max_cnt) {
                recv_task_id[i]++;
                if (recv_task_id[i] >= taskcnt) break;
                working_task++;
                cnt = 0;
                max_cnt = expected_recvcnt[ taskcnt * i + working_task ];
            }
        }

        recvcnt[i] = cnt;
    }

    #if LOG_LEVEL >= 3
    timer.stop_and_log("receive");
    #endif
}


void
prepare_supermer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int thr_per_task,
     int max_thr_membounded)
{

    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    Logger logger(comm);
    std::ostringstream rootlog;

    /* parallel settings */
    omp_set_nested(1);
    int ntasks = omp_get_max_threads() / thr_per_task;
    if (ntasks < 1) {
        ntasks = 1;
        thr_per_task = omp_get_max_threads();
    }
    /* for memory bounded operations in this stage, we use another number of threads */
    int nthr_membounded = std::min(omp_get_max_threads() , max_thr_membounded);

    size_t numreads = myreads.size();     /* Number of locally stored reads */
    // !! There may be a bug with large file and small rank num. in this case myread can read things wrong, resulting in 
    // !! Kmers with all A and C.


    #if LOG_LEVEL >= 2
    logger() << ntasks << " \t (thread per task: " << thr_per_task << ")";
    logger.flush("Task num:");
    logger() << nthr_membounded ;
    logger.flush("Thread count used for memory bounded operations:");
    #endif

    #if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
    #endif

    size_t readoffset = numreads;      
    MPI_Exscan(&numreads, &readoffset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
      

    /* data structure for storing the data of different threads */
    ParallelData data(nprocs, ntasks, nthr_membounded);

    /* find the destinations of each kmer */
    FindKmerDestinationsParallel(myreads, nthr_membounded, ntasks*nprocs, data);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer destination finding");
    timer.start();
    #endif

    /* log output for debugging */
    auto dest = data.get_my_destinations(0);

    for(int i = 0; i < dest.size(); i++) {
        for (int j = 0; j < dest[i].size(); j++) {
            logger() << dest[i][j] << " ";
        }
        logger() << std::endl;
    }
    logger.flush("Kmer destinations");

    /* encode the supermers */
    #pragma omp parallel num_threads(nthr_membounded)
    {
        int tid = omp_get_thread_num();
        auto& lengths = data.get_my_lengths(tid);
        auto& supermers = data.get_my_supermers(tid);
        auto& destinations = data.get_my_destinations(tid);
        auto& readids = data.get_my_readids(tid);

        const int MAX_SUPERMER_LEN = 256;

        SupermerEncoder encoder(lengths, supermers, MAX_SUPERMER_LEN);

        for (size_t i = 0; i < readids.size(); ++i) {
            encoder.encode(destinations[i], myreads[readids[i]]);
        }

    }

    /* log output for debugging */
    for (int tid = 0; tid < nthr_membounded; tid++) {
        std::vector<size_t> cur_locs;
        cur_locs.resize(nprocs * ntasks, 0);

        auto& lengths = data.get_my_lengths(tid);
        auto& supermers = data.get_my_supermers(tid);
        auto& destinations = data.get_my_destinations(tid);
        auto& readids = data.get_my_readids(tid);

        for(size_t i = 0; i < nprocs * ntasks; i++) {

            for(size_t j = 0; j < lengths[i].size(); j++) {
                logger() << "dst: "<< i << " length: " << lengths[i][j] << " ";
                logger() << "supermer: ";
                
                for(size_t k = 0; k < lengths[i][j]; k++) {
                    bool a = supermers[i][cur_locs[i] + k * 2];
                    bool b = supermers[i][cur_locs[i] + k * 2 + 1];
                    logger() << ((a && b) ? "T" : (a ?  "G" : (b ? "C" : "A")));
                    // 0 means A, 1 means C, 2 means G, and 3 means T.
                }   

                logger() << std::endl;

                cur_locs[i] += (lengths[i][j] * 2 + 6) % 8 * 8;
            }

        }
    }
    logger.flush("Supermer encoding");
    
    // return std::unique_ptr<KmerSeedBuckets>(recv_kmerseeds);
}


std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, int thr_per_task)
{

    Logger logger(MPI_COMM_WORLD);
    int myrank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    std::ostringstream rootlog;

    /* determine how many threads to use for each task */
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
    Timer timer(MPI_COMM_WORLD);
    timer.start();
    #endif

    /* sort the receiving vectors */
    #pragma omp parallel 
    {
        int tid = omp_get_thread_num();
        paradis::sort<KmerSeedStruct, TKmer::NBYTES>((*recv_kmerseeds)[tid].data(), (*recv_kmerseeds)[tid].data() + task_seedcnt[tid], thr_per_task);
    

#if LOG_LEVEL >= 3
    }
    timer.stop_and_log("Shared memory parallel K-mer sorting");
    print_mem_log(nprocs, myrank, "After Sorting");
    timer.start();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
#endif

        /* do a linear scan and filter out the kmers we want */
        kmerlists[tid].reserve(uint64_t(task_seedcnt[tid] / LOWER_KMER_FREQ));
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
                    int& count              = std::get<1>(entry);

                    count = cur_kmer_cnt;
                    kmer = last_mer;

                    valid_kmer[tid]++;
                }

                cur_kmer_cnt = 1;
                last_mer = cur_mer;
            }
        }

        /* deal with the last kmer of a task */
        if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {
            kmerlists[tid].push_back(KmerListEntry());
            KmerListEntry& entry         = kmerlists[tid].back();

            TKmer& kmer             = std::get<0>(entry);
            int& count              = std::get<1>(entry);

            count = cur_kmer_cnt;
            kmer = last_mer;


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
    print_mem_log(nprocs, myrank, "After Counting");
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

// unsigned short should be enough here
void FindKmerDestinationsParallel(const DnaBuffer& myreads, int nthreads, int tot_tasks, ParallelData& data) {
    
    assert(nthreads > 0);

    // Hope the static is what we want 
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (size_t i = 0; i < myreads.size(); ++i)
    {
        int tid = omp_get_thread_num();
        std::cout<<"tid: "<<tid<<" is responsible for read "<<i<<"\n";

        auto &dest = data.register_new_destination(tid, i);
        if (myreads[i].size() < KMER_SIZE)
            continue;
        dest.reserve(myreads[i].size() - KMER_SIZE + 1);

        std::vector<TMmer> repmers = TMmer::GetRepMmers(myreads[i]);

        Minimizer_Deque deque;
        int head_pos = 0;

        // Initialize the deque
        for(; head_pos < KMER_SIZE - MINIMIZER_SIZE; head_pos++) {
            deque.insert_minimizer(repmers[head_pos].GetHash(), head_pos);
        }
        int tail_pos = head_pos - KMER_SIZE + MINIMIZER_SIZE - 1;

        // Start the main loop
        for(; head_pos < repmers.size(); head_pos++, tail_pos++) {
            deque.insert_minimizer(repmers[head_pos].GetHash(), head_pos);
            // std::cout<<"Hash is"<<repmers[head_pos].GetHash()<<" "<<deque.get_current_minimizer()<<std::endl;
            deque.remove_minimizer(tail_pos);
            dest.push_back(GetMinimizerOwner(deque.get_current_minimizer(), tot_tasks));
        }

        // Need To Check The Code Above!
    }
}


inline int GetMinimizerOwner(const uint64_t& hash, int tot_tasks) {
    // Need to check if this gives equal distribution
    return static_cast<int>(hash % tot_tasks);
}


inline int GetKmerOwner(const TKmer& kmer, int nprocs) {
    uint64_t myhash = kmer.GetHash();
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / std::numeric_limits<uint64_t>::max();
    assert(owner >= 0 && owner < static_cast<int>(nprocs));
    return static_cast<int>(owner);
}


void print_kmer_histogram(const KmerList& kmerlist, MPI_Comm comm) {
    #if LOG_LEVEL >= 2

    Logger logger(comm);
    int maxcount = std::accumulate(kmerlist.cbegin(), kmerlist.cend(), 0, [](int cur, const auto& entry) { return std::max(cur, std::get<1>(entry)); });

    MPI_Allreduce(MPI_IN_PLACE, &maxcount, 1, MPI_INT, MPI_MAX, comm);

    std::vector<int> histo(maxcount+1, 0);

    for(size_t i = 0; i < kmerlist.size(); ++i)
    {
        int cnt = std::get<1>(kmerlist[i]);
        assert(cnt >= 1);
        histo[cnt]++;
    }

    MPI_Allreduce(MPI_IN_PLACE, histo.data(), maxcount+1, MPI_INT, MPI_SUM, comm);

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        std::cout << "#count\tnumkmers" << std::endl;

        for (int i = 1; i < histo.size(); ++i)
        {
            if (histo[i] > 0)
            {
                std::cout << i << "\t" << histo[i] << std::endl;
            }
        }
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
    #endif
}
