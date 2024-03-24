#include "kmerops.hpp"
#include "logger.hpp"
#include "dnaseq.hpp"
#include "timer.hpp"
#include "paradissort.hpp"
#include "memcheck.hpp"
#include "supermer.hpp"
#include "compiletime.h"
#include "raduls.h"
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
     int thr_per_worker,
     int max_thr_membounded)
{

    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    Logger logger(comm);

    /* parallel settings */
    omp_set_nested(1);
    int ntasks = omp_get_max_threads() / thr_per_worker * AVG_TASK_PER_WORKER;
    if (ntasks < 1) {
        ntasks = AVG_TASK_PER_WORKER;
        thr_per_worker = omp_get_max_threads();
    }
    /* for memory bounded operations in this stage, we use another number of threads */
    int nthr_membounded = std::min(omp_get_max_threads() , max_thr_membounded);

    size_t numreads = myreads.size();     /* Number of locally stored reads */
    // !! There may be a bug with large file and small rank num. in this case myread can read things wrong, resulting in 
    // !! Kmers with all A and C.


#if LOG_LEVEL >= 2
    logger() << ntasks << " \t (thread per worker: " << thr_per_worker << ")";
    logger.flush("Task num:", 0);
    logger() << nthr_membounded ;
    logger.flush("Thread count used for memory bounded operations:", 0);
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
    timer.stop_and_log("(Inc) K-mer destination finding");
    timer.start();
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

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Supermer encoding");
#endif

    return data;
}



std::unique_ptr<KmerSeedBuckets> exchange_supermer(ParallelData& data, MPI_Comm comm, TaskDispatcher& dispatcher)
{
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
    Logger logger(comm);

#if LOG_LEVEL >= 3
Timer timer(comm);
timer.start();
#endif


    auto local_tasksz = data.get_local_tasksz();
    dispatcher.balanced_dispatch(comm, local_tasksz);
    // dispatcher.plain_dispatch();
    int mytasks = dispatcher.get_taskid(myrank).size();

#if LOG_LEVEL >= 3
timer.stop_and_log("(Inc) Task dispatch");
timer.start();
#endif

    std::vector<size_t> send_counts(nprocs * data.ntasks, 0);
    std::vector<size_t> recv_counts(nprocs * mytasks, 0);

    std::vector<int32_t> scounts(nprocs, 0);
    std::vector<int32_t> sdispls(nprocs, 0);
    std::vector<int32_t> rcounts(nprocs, 0);
    std::vector<int32_t> rdispls(nprocs, 0);


    uint32_t offset = 0;
    for (int i = 0; i < nprocs; i++) {
        auto taskids = dispatcher.get_taskid(i);
        scounts[i] = taskids.size();
        sdispls[i] = offset;
        for (int j = 0; j < taskids.size(); j++) {
            send_counts[offset + j] = data.get_supermer_cnt(taskids[j]);
        }
        offset += taskids.size();
    }


    for (int i = 0; i < nprocs; i++) {
        rdispls[i] = i * mytasks;
        rcounts[i] = mytasks;
    }


    MPI_Alltoallv(send_counts.data(), scounts.data(), sdispls.data(), MPI_UNSIGNED_LONG_LONG, 
            recv_counts.data(), rcounts.data(), rdispls.data(), MPI_UNSIGNED_LONG_LONG, comm);

    std::vector<std::vector<uint32_t>> length;

#if LOG_LEVEL >= 3
timer.stop_and_log("(Inc) Task Size information exchange");
timer.start();
#endif

    LengthExchanger length_exchanger(comm, MAX_SEND_BATCH, sizeof(uint32_t), data.nthr_membounded, data.lengths, recv_counts, length, dispatcher);
    length_exchanger.initialize();

    while(length_exchanger.status != LengthExchanger::Status::BATCH_DONE) {
        length_exchanger.progress();
    }

#if LOG_LEVEL >= 3
timer.stop_and_log("(Inc) Length exchange");
length_exchanger.print_stats();
timer.start();
#endif

    KmerSeedBuckets* bucket = new KmerSeedBuckets(mytasks);
    SupermerExchanger supermer_exchanger(comm, MAX_SEND_BATCH, MAX_SUPERMER_LEN, data.nthr_membounded, data.lengths, data.supermers, recv_counts, length, *bucket, dispatcher); 
    
    supermer_exchanger.initialize();
    while (supermer_exchanger.status != SupermerExchanger::Status::BATCH_DONE)
    {
        supermer_exchanger.progress();
    }

#if LOG_LEVEL >= 3
timer.stop_and_log("(Inc) Supermer exchange");
supermer_exchanger.print_stats();
#endif
    

    return std::unique_ptr<KmerSeedBuckets>(bucket);
}

int get_ppn() {
    char* ppnc = std::getenv("SLURM_TASKS_PER_NODE");
    if (ppnc == NULL) {
        return -1;
    }
    // read until meet the first non-digit character
    int ppn = 0;
    while (*ppnc != '\0' && *ppnc >= '0' && *ppnc <= '9') {
        ppn = ppn * 10 + *ppnc - '0';
        ppnc++;
    }
    return ppn;
}


std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, TaskDispatcher& dispatcher, int thr_per_worker)
{

    Logger logger(MPI_COMM_WORLD);
    int myrank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    std::ostringstream rootlog;

    /* determine how many threads to use for each task */
    int nworkers = omp_get_max_threads() / thr_per_worker;
    if (nworkers < 1){
        nworkers = 1;
        thr_per_worker = omp_get_max_threads();
    }

    auto tasks = dispatcher.get_taskid(myrank);
    int mytasks = tasks.size();

    assert(mytasks == recv_kmerseeds->size());

    omp_set_nested(1);
    omp_set_num_threads(nworkers);

    uint64_t task_seedcnt[mytasks];
    uint64_t task_seedtot = 0;
    uint64_t valid_kmer[mytasks];
    KmerList kmerlists[mytasks];


    for (int i = 0; i < mytasks; i++) {
        task_seedcnt[i] = (*recv_kmerseeds)[i].size();
        task_seedtot += task_seedcnt[i];
    }

    
#if LOG_LEVEL >= 3
    logger() <<"worker num: "<< nworkers << " \tthread per worker: " << thr_per_worker << std::endl;
    logger() <<"task num: "<<mytasks<<" \ttask size: ";
    for (int i = 0; i < mytasks; i++) {
        logger() << task_seedcnt[i] << " ";
    }
    logger.flush("Parallel tasks for kmer sorting and counting:");

    size_t mytotseed = std::accumulate(task_seedcnt, task_seedcnt + mytasks, (uint64_t)0);
    uint64_t all_seedcnt[nprocs];
    MPI_Gather(&mytotseed, 1, MPI_UNSIGNED_LONG_LONG, all_seedcnt, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    // rank 0 calculate the standard deviation
    if (myrank==0) {
        uint64_t sum = 0;
        for (int i = 0; i < nprocs; i++) {
                sum += all_seedcnt[i];
        }
        uint64_t mean = (sum) / (nprocs);
        uint64_t stddev = 0;
        for (int i = 0; i < nprocs; i++) {
            stddev += (all_seedcnt[i] - mean) * (all_seedcnt[i] - mean);
        }
        stddev = sqrt(stddev / (nprocs));
        logger() << "mean: " << mean << " \t\t std_dev: " << stddev << " \t\t ratio: " << (double)stddev / mean;
    }
    logger.flush("Task size count distribution statistics:", 0);

    Timer timer(MPI_COMM_WORLD);
    timer.start();
#endif


    /* choose the sort algorithm */ 

    int sort = 0;

    if( SORT == 1 ) {
        sort = 1;
        logger() << "Using PARADIS for sorting.";
    } else if( SORT == 2 ) {
        sort = 2;
        logger() << "Using RADULS for sorting.";
    } else {
        /* SORT == 0, decide upon available memory */
        int ppn = get_ppn(); 
        size_t memfree_kb; 
        if (get_free_memory_kb(&memfree_kb) == 1 || ppn == -1) {
            logger() << "Warning: Could not get free memory or Process per node. Default to PARADIS.";
            sort = 1;
        } else {
            size_t memfree = memfree_kb * 921 / ppn ;   // 1024 * 0.9 = 921
            if (memfree >= task_seedtot * TKmer::NBYTES) {
                sort = 2;
                logger() << "Enough memory available. Using RADULS for sorting. Average available mem (Bytes): "<<memfree<<"  Kmer Seeds:"<<task_seedtot<<"  Processes Per Node:"<<ppn<< " Bytes per Kmer:" << TKmer::NBYTES;
            } else {
                logger() << "Not enough memory available. Using PARADIS for sorting.";
                sort = 1;
            }
        }
    }
    logger.flush("Sort algorithm decision:");
    print_mem_log(nprocs, myrank, "Before Sorting");

    size_t start_pos[mytasks];

    double times[nworkers];
    int unfinished_workers = std::min(nworkers, mytasks);
    int total_threadnum = thr_per_worker * nworkers;
    bool tasks_processed[mytasks];
    for (int i = 0; i < mytasks; i++) {
        tasks_processed[i] = false;
    }

    /* sort the receiving vectors */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        // int current_task_idx = tid;
        TimerLocal tm;
        tm.start();

        while(true) {
            int current_task = -1;
            bool finished = false;
            #pragma omp critical
            {
                for (int i = 0; i < mytasks; i++) {
                    if (!tasks_processed[i]) {
                        current_task = i;
                        tasks_processed[i] = true;
                        break;
                    }
                }
                if (current_task == -1) {
                    finished = true;
                }
            }
            if (finished) {
                break;
            }

            int thr_per_worker = total_threadnum / unfinished_workers;

            #pragma omp critical
            {
                logger()<< "Task "<<current_task<<" Sorting with "<<thr_per_worker<<" threads on worker "<<tid<<std::endl;
            }
            if (sort == 1){
                start_pos[current_task] = 0;
                paradis::sort<KmerSeedStruct, TKmer::NBYTES>((*recv_kmerseeds)[current_task].data(), (*recv_kmerseeds)[current_task].data() + task_seedcnt[current_task], thr_per_worker);
            } else {
                uint8_t* tmp_arr = new uint8_t[task_seedcnt[current_task] * TKmer::NBYTES + 256];
                uint8_t* tmp = tmp_arr + 256 - (size_t)tmp_arr % 256;
                raduls::CleanTmpArray(tmp, task_seedcnt[current_task], TKmer::NBYTES, thr_per_worker);

                // RADULS needs padding. reserved in bucket assistant
                uint8_t* start = (uint8_t*)(*recv_kmerseeds)[current_task].data();
                int cnt = 0;
                while( (size_t)start % 256 != 0) {
                    start += TKmer::NBYTES;
                    (*recv_kmerseeds)[current_task].push_back((*recv_kmerseeds)[current_task][cnt]);
                    cnt++;
                }
                start_pos[current_task] = cnt;

                raduls::RadixSortMSD(start, tmp, task_seedcnt[current_task], TKmer::NBYTES, TKmer::NBYTES, thr_per_worker);
                delete[] (tmp_arr);
            }

        }
        times[tid] = tm.stop();
        #pragma omp critical
        {
            unfinished_workers--;
        }
    }

    logger.flush("Sorting time for tasks:");


    
#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Shared memory parallel K-mer sorting");
    print_mem_log(nprocs, myrank, "After Sorting");
    logger() << "Sorting time: ";
    for (int i = 0; i < nworkers; i++) {
        logger() << times[i] << " ";
    }
    logger.flush("Sorting time for tasks:");

    timer.start();
#endif

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int current_task_idx = tid;

        while(current_task_idx < mytasks) {
            int current_task = current_task_idx;

            /* do a linear scan and filter out the kmers we want */
            
            kmerlists[current_task].reserve(uint64_t(task_seedcnt[current_task] / LOWER_KMER_FREQ));
            TKmer last_mer = (*recv_kmerseeds)[current_task][start_pos[current_task]].kmer;
            uint64_t cur_kmer_cnt = 1;
            valid_kmer[current_task] = 0;

            for(size_t idx = start_pos[current_task] + 1; idx < task_seedcnt[current_task] + start_pos[current_task]; idx++) {
                TKmer cur_mer = (*recv_kmerseeds)[current_task][idx].kmer;
                if (cur_mer == last_mer) {
                    cur_kmer_cnt++;
                } else {
                    // if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ)
                    if (cur_kmer_cnt >= LOWER_KMER_FREQ ) {

                        kmerlists[current_task].push_back(KmerListEntry());
                        KmerListEntry& entry    = kmerlists[current_task].back();

                        TKmer& kmer             = std::get<0>(entry);
                        int& count              = std::get<1>(entry);

                        count = cur_kmer_cnt;
                        kmer = last_mer;

                        valid_kmer[current_task]++;
                    }

                    cur_kmer_cnt = 1;
                    last_mer = cur_mer;
                }
            }

            /* deal with the last kmer of a task */
            // if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) 
            if (cur_kmer_cnt >= LOWER_KMER_FREQ ) {
                kmerlists[current_task].push_back(KmerListEntry());
                KmerListEntry& entry         = kmerlists[current_task].back();

                TKmer& kmer             = std::get<0>(entry);
                int& count              = std::get<1>(entry);

                count = cur_kmer_cnt;
                kmer = last_mer;


                valid_kmer[current_task]++;
            }
            current_task_idx += nworkers;
        }
    }

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) K-mer counting");
    timer.start();
#endif

#if LOG_LEVEL >= 2
    for (int i = 0; i < mytasks; i++) {
        logger() << valid_kmer[i] << " ";
    }
    logger.flush("Valid kmer from tasks:");
    print_mem_log(nprocs, myrank, "After Counting");
#endif


    uint64_t valid_kmer_total = std::accumulate(valid_kmer, valid_kmer + mytasks, (uint64_t)0);
    KmerList* kmerlist = new KmerList();
    kmerlist->reserve(valid_kmer_total);

    // yfli: actually we don't need to copy it into one single kmerlist
    for (int i = 0; i < mytasks; i++) {
        kmerlist->insert(kmerlist->end(), kmerlists[i].begin(), kmerlists[i].end());
    }

    logger() << valid_kmer_total;
    logger.flush("Valid kmer for process:");

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) K-mer copying");
#endif

#if LOG_LEVEL >= 4
    // write kmer list to a file
    std::ostringstream kmerfilename;
    kmerfilename << "/pscratch/sd/y/yfli03/log/kmerlist_" << myrank << ".txt";
    
    std::ofstream kmerfile;
    kmerfile.open(kmerfilename.str());
    for (size_t i = 0; i < kmerlist->size(); i++) {
        kmerfile << std::get<0>((*kmerlist)[i]) << " " << std::get<1>((*kmerlist)[i]) << std::endl;
    }

    kmerfile.close();
#endif

    return std::unique_ptr<KmerList>(kmerlist);
}


void FindKmerDestinationsParallel(const DnaBuffer& myreads, int nthreads, int tot_tasks, ParallelData& data) {
    
    assert(nthreads > 0);

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

        /* initialize the deque */
        for(; head_pos < KMER_SIZE - MINIMIZER_SIZE; head_pos++) {
            deque.insert_minimizer(repmers[head_pos].GetHash(), head_pos);
        }
        int tail_pos = head_pos - KMER_SIZE + MINIMIZER_SIZE - 1;

        /* start the main loop */
        for(; head_pos < repmers.size(); head_pos++, tail_pos++) {
            deque.insert_minimizer(repmers[head_pos].GetHash(), head_pos);
            deque.remove_minimizer(tail_pos);
            dest.push_back(GetMinimizerOwner(deque.get_current_minimizer(), tot_tasks));
        }

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
