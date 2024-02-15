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


ParallelData
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
#if LOG_LEVEL >= 4
    auto dest = data.get_my_destinations(0);

    for(int i = 0; i < dest.size(); i++) {
        for (int j = 0; j < dest[i].size(); j++) {
            logger() << dest[i][j] << " ";
        }
        logger() << std::endl;
    }
    logger.flush("Kmer destinations");
#endif

    /* encode the supermers */
    #pragma omp parallel num_threads(nthr_membounded)
    {
        int tid = omp_get_thread_num();
        auto& lengths = data.get_my_lengths(tid);
        auto& supermers = data.get_my_supermers(tid);
        auto& destinations = data.get_my_destinations(tid);
        auto& readids = data.get_my_readids(tid);

        SupermerEncoder encoder(lengths, supermers, MAX_SUPERMER_LEN);

        for (size_t i = 0; i < readids.size(); ++i) {
            encoder.encode(destinations[i], myreads[readids[i]]);
        }

    }

#if LOG_LEVEL >= 4
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
                    bool a = supermers[i][cur_locs[i] + k / 4] & (1 << (7 - (k%4)*2));
                    bool b = supermers[i][cur_locs[i] + k / 4] & (1 << (6 - (k%4)*2));
                    logger() << ((a && b) ? "T" : (a ?  "G" : (b ? "C" : "A")));
                    /* 0 means A, 1 means C, 2 means G, and 3 means T. */
                }   

                logger() << std::endl;

                cur_locs[i] += cnt_bytes(lengths[i][j]);
            }

        }
    }
    logger.flush("Supermer encoding");
#endif

    return data;
}



std::unique_ptr<KmerSeedBuckets> exchange_supermer(ParallelData& data, MPI_Comm comm)
{
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
    Logger logger(comm);

    int ntasks = data.ntasks;
    // Calculate the count of supermers to be sent to each task
    std::vector<size_t> send_counts(nprocs * ntasks, 0);
    std::vector<size_t> recv_counts(nprocs * ntasks, 0);

    for (int i = 0; i < nprocs; i++) {
        for (int j = 0; j < ntasks; j++) 
            send_counts[i * ntasks + j] = data.get_supermer_cnt(i, j);
    }

    MPI_Alltoall(send_counts.data(), ntasks, MPI_UNSIGNED_LONG_LONG, recv_counts.data(), ntasks, MPI_UNSIGNED_LONG_LONG, comm);

    std::vector<std::vector<uint32_t>> length;

#if LOG_LEVEL >= 3
Timer timer(comm);
timer.start();
#endif

    // Exchange the lengths
    LengthExchanger length_exchanger(comm, ntasks, MAX_SEND_BATCH, sizeof(uint32_t), data.nthr_membounded, data.lengths, recv_counts, length);
    length_exchanger.initialize();

    while(length_exchanger.status != LengthExchanger::Status::BATCH_DONE) {
        length_exchanger.progress();
    }
#if LOG_LEVEL >= 3
timer.stop_and_log("Length exchange");
length_exchanger.print_stats();
timer.start();
#endif

#if LOG_LEVEL >= 4
    for (int i = 0; i < length.size(); i++) {
        logger() << "recvbuf " << i << std::endl;
        for (int j = 0; j < length[i].size(); j++) {
            logger() << length[i][j] << " ";
        }
        logger() << std::endl;
    }
    
    logger.flush("Received lengths");
#endif

    KmerSeedBuckets* bucket = new KmerSeedBuckets(ntasks);
    // Exchange the supermers
    SupermerExchanger supermer_exchanger(comm, ntasks, MAX_SEND_BATCH, MAX_SUPERMER_LEN, data.nthr_membounded, data.lengths, data.supermers, recv_counts, length, *bucket); 
    
    supermer_exchanger.initialize();
    while (supermer_exchanger.status != SupermerExchanger::Status::BATCH_DONE)
    {
        supermer_exchanger.progress();
    }

#if LOG_LEVEL >= 3
timer.stop_and_log("Supermer exchange");
supermer_exchanger.print_stats();
#endif
    

    return std::unique_ptr<KmerSeedBuckets>(bucket);
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
