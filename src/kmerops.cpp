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
#include <stdexcept>



std::shared_ptr<TaskManager>
prepare_supermer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int avg_thr_per_worker,
     int max_thr_m)
{

    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    Logger logger(comm);

    /* parallel settings */
    omp_set_nested(1);
    // The '-1' is an optimization for heavy hitter strategy
    int avg_tasks = omp_get_max_threads() / avg_thr_per_worker * AVG_TASK_PER_WORKER - 1;
    if (avg_tasks < 1) {
        avg_tasks = AVG_TASK_PER_WORKER;
    }

    /* for memory bounded operations in this stage, we use another number of threads */
    int nthr_m = std::min(omp_get_max_threads() , max_thr_m);

    size_t numreads = myreads.size();     /* Number of locally stored reads */
    // !! There may be a bug with large file and small rank num. in this case myread can read things wrong, resulting in 
    // !! Kmers with all A and C.


#if LOG_LEVEL >= 2
    logger() << "Total task count: " << avg_tasks * nprocs << std::endl
             << "  Total thread count: " << omp_get_max_threads() << std::endl
             << "  Total thread count for memory bounded operations: " << nthr_m;
    logger.flush("Initialize info:", 0);
#endif

#if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
#endif

    size_t readoffset = numreads;      
    MPI_Exscan(&numreads, &readoffset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
      
    /* data structure for storing the data of different threads */
    ParallelData data(nthr_m);

    /* find the destinations of each kmer */
    FindKmerDestinationsParallel(myreads, nthr_m, avg_tasks * nprocs, data);

    /* initialize the tasks */
    std::vector<std::shared_ptr<ScatteredTask>> sca_tasks;
    for(int i = 0; i < avg_tasks * nprocs; i++) {
        sca_tasks.push_back(std::make_shared<ScatteredSupermers>(i, nthr_m));
    }


#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) K-mer destination finding");
    timer.start();
#endif

    /* encode the supermers */
    #pragma omp parallel num_threads(nthr_m)
    {
        int tid = omp_get_thread_num();
        SupermerEncoder encoder(sca_tasks, tid, MAX_SUPERMER_LEN);

        auto& destinations = data.get_my_destinations(tid);
        auto& readids = data.get_my_readids(tid);

        for (size_t i = 0; i < readids.size(); ++i) {
            encoder.encode(destinations[i], myreads[readids[i]]);
        }

    }

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Supermer encoding");
#endif

    HeavyHitterClassifier* classifier = new HeavyHitterClassifier();
    BalancedDispatcher* dispatcher = new BalancedDispatcher();


    auto tm = std::make_shared<TaskManager>(classifier, dispatcher, comm, MAX_SEND_BATCH, MAX_SEND_BATCH - std::max((uint64_t)MAX_SUPERMER_LEN / 4, sizeof(KmerSeedStruct)), nthr_m);
    tm->init_scattered_tasks(sca_tasks);

    return tm;
}



void exchange_supermer(std::shared_ptr<TaskManager> task_manager, MPI_Comm comm, int avg_thr_per_worker)
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

    if (avg_thr_per_worker > omp_get_max_threads()) {
        avg_thr_per_worker = omp_get_max_threads();
    }

    TaskManager& tm = *task_manager;

    tm.classify_tasks();

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task classification");
    timer.start();
#endif

    tm.preprocess_tasks(avg_thr_per_worker);

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task preprocessing");
    timer.start();
#endif

    tm.dispatch();

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task dispatch");
    timer.start();
#endif

    // Create my gathered tasks
    tm.init_exchange_groups();  

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task exchange group initialization");   
    timer.start();
#endif  

    // Start exchange stage 1
    tm.exchange(1);

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task exchange stage 1");
    timer.start();
#endif

    // Start exchange stage 2
    tm.exchange(2);

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) Task exchange stage 2");
#endif

    return;

}


std::unique_ptr<KmerListS>
filter_kmer(std::shared_ptr<TaskManager> task_manager, MPI_Comm comm, int avg_thr_per_worker)
{
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    Logger logger(comm);

#if LOG_LEVEL >= 3
    Timer timer(comm);
    timer.start();
    TimerLocal tl;
    tl.start();
#endif


    /* determine how many threads to use for each task */
    int nworkers = omp_get_max_threads() / avg_thr_per_worker;
    if (nworkers < 1){
        nworkers = 1;
        avg_thr_per_worker = omp_get_max_threads();
    }

    omp_set_nested(1);

    TaskManager& tm = *task_manager;

    size_t gathered_task_bytes = tm.get_gathered_task_bytes();

    /* choose the sort algorithm */ 
    int sort = sort_decision(gathered_task_bytes, logger);
    logger.flush("Sort algorithm decision:");

    tm.process_tasks(nworkers, avg_thr_per_worker, sort);

#if LOG_LEVEL >= 3
    logger() << tl.stop() << std::endl;
    timer.stop_and_log("(Inc) K-mer counting");
    logger.flush("K-mer counting details:");
    timer.start();
#endif

    KmerListS* kmerlist = new KmerListS();
    tm.copy_results(*kmerlist);

#if LOG_LEVEL >= 3
    timer.stop_and_log("(Inc) K-mer copying");
#endif

    return std::unique_ptr<KmerListS>(kmerlist);
}


ScatteredSupermers::ScatteredSupermers(int taskid, int nthr) : 
    ScatteredTask(taskid, 1), nthr_(nthr)
    {
        supermers_.resize(nthr_);
        lengths_.resize(nthr_);
    }

size_t ScatteredSupermers::get_size_kmer() {
    size_t size = 0;
    for (int i = 0; i < nthr_; i++) {
        size += std::accumulate(lengths_[i].begin(), lengths_[i].end(), 0);
        size -= ( KMER_SIZE - 1 ) * lengths_[i].size() ;
    }
    return size;
}

size_t ScatteredSupermers::get_size()  {
    size_t size = 0;
    for (int i = 0; i < nthr_; i++) {
        size += lengths_[i].size();
    }
    return size;
}

size_t ScatteredSupermers::get_size_bytes()  {
    size_t size = 0;
    for (int i = 0; i < nthr_; i++) {
        size += lengths_[i].size() * sizeof(uint32_t);
    }
    for (int i = 0; i < nthr_; i++) {
        size += supermers_[i].size();
    }
    return size;
}

void ScatteredSupermers::init_exchange(int stage)  {
    working_thr_ = 0;
    working_idx_ = 0;
    working_buf_idx_ = 0;
}

bool ScatteredSupermers::write_to_buffer_stage1(uint8_t* buffer, size_t &current_count, size_t max_count) {
    while (working_thr_ < nthr_ && current_count <= max_count) {

        auto& lengths = lengths_[working_thr_];
        size_t copy_num = std::min((max_count - current_count) / LEN_SIZE + 1, (lengths.size() - working_idx_));
        memcpy(buffer + current_count, lengths.data() + working_idx_, copy_num * LEN_SIZE);
        current_count += copy_num * LEN_SIZE;

        working_idx_ += copy_num;
        if (working_idx_ == lengths.size()) {
            working_idx_ = 0;
            working_thr_++;
        }

    }
    return working_thr_ == nthr_;
}

bool ScatteredSupermers::write_to_buffer_stage2(uint8_t* buffer, size_t &current_count, size_t max_count) {
    while(working_thr_ < nthr_ && current_count <= max_count) {

        auto& supermers = supermers_[working_thr_];
        auto& lengths = lengths_[working_thr_];

        size_t copy_len = 0;
        size_t idx = working_idx_;

        while(current_count + copy_len <= max_count && idx < lengths.size()) {
            copy_len += cnt_bytes(lengths[idx]);
            idx++;
        }

        memcpy(buffer + current_count, supermers.data() + working_buf_idx_, copy_len);
        current_count += copy_len;

        working_idx_ = idx;
        working_buf_idx_ += copy_len;
        if (working_idx_ == lengths.size()) {
            working_idx_ = 0;
            working_buf_idx_ = 0;
            working_thr_++;
        }

    }
    return working_thr_ == nthr_;
}

bool ScatteredSupermers::write_to_buffer(uint8_t* buffer, size_t& current_count, size_t max_count, int stage)  {
    if (stage == 1) {
        return write_to_buffer_stage1(buffer, current_count, max_count);
    } else if (stage == 2) {
        return write_to_buffer_stage2(buffer, current_count, max_count);
    }
    throw std::runtime_error("Unknown stage.");
}




ScatteredKmerList::ScatteredKmerList(std::shared_ptr<ScatteredSupermers> supermer_task, int thr_per_worker, size_t stage1_size) : 
        working_idx_(0), stage1_size_(stage1_size),  ScatteredTask(supermer_task->get_taskid(), 2) {

    std::vector<KmerSeedStruct> kmerseeds;

    int nthr = supermer_task->get_nthr();
    size_t total_len = supermer_task->get_size_kmer();
    kmerseeds.reserve(total_len + MAX_SUPERMER_LEN / TKmer::NBYTES);

    /* extract all the kmers from supermers */
    for (int i = 0; i < nthr; i++) {
        size_t idx = 0;
        auto& lengths = supermer_task->get_length_buffer(i);
        auto& supermers = supermer_task->get_supermer_buffer(i);

        for (size_t j = 0; j < lengths.size(); j++) {
            size_t len = lengths[j];
            size_t len_bytes = cnt_bytes(len);
            auto seq = DnaSeq(len, supermers.data() + idx);
            auto repmers = TKmer::GetRepKmers(seq);
            for (int k = 0; k < len - KMER_SIZE + 1; k++) {
                kmerseeds.emplace_back(repmers[k]);
            }
            idx += len_bytes;
        }
    }

    assert(kmerseeds.size() == total_len);
    size_t start_pos;
    sort_task(kmerseeds, 2, thr_per_worker, start_pos, total_len);
    size_t valid_kmer;
    count_sorted_kmers(kmerseeds, kmerlist_, start_pos, total_len, valid_kmer, false);

    stage1_size_ = std::min(stage1_size_, kmerlist_.size());
}


bool ScatteredKmerList::write_to_buffer(uint8_t* buffer, size_t& current_count, size_t max_count, int stage)  {
    size_t end_idx = kmerlist_.size();
    if (stage == 1) {
        end_idx = stage1_size_;
    }

    size_t copy_num = std::min((max_count - current_count) / sizeof(KmerListEntryS) + 1, end_idx - working_idx_);
    memcpy(buffer + current_count, kmerlist_.data() + working_idx_, copy_num * sizeof(KmerListEntryS));
    current_count += copy_num * sizeof(KmerListEntryS);
    working_idx_ += copy_num;

    return working_idx_ == end_idx;
    
}


GatheredTask::GatheredTask(int taskid, int tasktype, std::vector<size_t> count) : 
        taskid_(taskid), tasktype_(tasktype), count_(count)
{
    working_idx_.resize(count.size(), 0);
    start_idx_.resize(count.size() + 1, 0);
    for (int i = 1; i <= count.size(); i++) {
        start_idx_[i] = start_idx_[i-1] + count[i-1];
        working_idx_[i-1] = start_idx_[i-1];
    }

}


void GatheredSupermer::init_exchange(int stage)  {
    if (stage == 1){
        length_.resize(std::accumulate(count_.begin(), count_.end(), 0));
        return;
    }

    if (stage == 2) {
        kmer_count_.resize(count_.size(), 0);
        for (int i = 0; i < count_.size(); i++) {
            kmer_count_[i] = std::accumulate(length_.begin() + start_idx_[i], length_.begin() + start_idx_[i+1], 0);
            kmer_count_[i] -= (KMER_SIZE - 1) * count_[i];

        }

        working_kmer_idx_.resize(kmer_count_.size(), 0);
        working_idx_[0] = 0;

        for (int i = 1; i < kmer_count_.size(); i++) {
            working_idx_[i] = start_idx_[i];
            working_kmer_idx_[i] = working_kmer_idx_[i-1] + kmer_count_[i-1];
        }

        size_t other_size = 0;
        for(int i = 0; i < count_.size(); i++) {
            other_size += kmer_count_[i];
        }

        kmer_.reserve(std::accumulate(kmer_count_.begin(), kmer_count_.end(), 0) + PAD_SIZE / sizeof(KmerSeedStruct) + 1);
        kmer_.resize(std::accumulate(kmer_count_.begin(), kmer_count_.end(), 0));

        return;
    }
    throw std::runtime_error("Unknown stage.");
}

bool GatheredSupermer::receive_from_buffer_stage1(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count) {
    size_t copy_num = std::min((max_count - current_count) / LEN_SIZE + 1, start_idx_[src_process+1] - working_idx_[src_process]);
    memcpy(length_.data() + working_idx_[src_process], buffer + current_count, copy_num * LEN_SIZE);
    current_count += copy_num * LEN_SIZE;
    working_idx_[src_process] += copy_num;

    return working_idx_[src_process] == start_idx_[src_process+1];
}

bool GatheredSupermer::receive_from_buffer_stage2(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count) {
    size_t current_kmer_idx = working_kmer_idx_[src_process];
    size_t current_idx = working_idx_[src_process];

    size_t ccnt = current_count;

    while(ccnt <= max_count && current_idx < start_idx_[src_process+1]) {

        auto seq = DnaSeq(length_[current_idx], buffer + ccnt);
        auto repmers = TKmer::GetRepKmers(seq);

        for (int i = 0; i < length_[current_idx] - KMER_SIZE + 1; i++) {
            kmer_[current_kmer_idx++] = KmerSeedStruct(repmers[i]);
        }

        ccnt += cnt_bytes(length_[current_idx]);
        current_idx++;
    }

    working_idx_[src_process] = current_idx;
    working_kmer_idx_[src_process] = current_kmer_idx;
    current_count = ccnt;


    return current_idx == start_idx_[src_process+1];
}


bool GatheredSupermer::receive_from_buffer(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count, int stage)  {

    if (stage == 1) {
        return receive_from_buffer_stage1(buffer, src_process, current_count, max_count);
    } else if (stage == 2) {
        return receive_from_buffer_stage2(buffer, src_process, current_count, max_count);
    }
    
    throw std::runtime_error("Unknown stage.");
}

void GatheredSupermer::process(int nthr, int sort)  {
    size_t start_pos;
    size_t total_kmer = std::accumulate(kmer_count_.begin(), kmer_count_.end(), 0);
    sort_task(kmer_, sort, nthr, start_pos, total_kmer);
    size_t valid_kmer;
    count_sorted_kmers(kmer_, filtered_kmerlist_, start_pos, total_kmer, valid_kmer, true);
}




GatheredKmerList::GatheredKmerList(int taskid, int tasktype, std::vector<size_t> count, size_t stage1_size) : 
    GatheredTask(taskid, tasktype, count)
{
    stage1_size_.resize(count.size(), stage1_size);
    for (int i = 0; i < count.size(); i++) {
        stage1_size_[i] = std::min(stage1_size, count[i]);
    }

    kmerlist_.reserve(std::accumulate(count.begin(), count.end(), 0) + PAD_SIZE / sizeof(KmerListEntryS) + 1);
    kmerlist_.resize(std::accumulate(count.begin(), count.end(), 0));
}

bool GatheredKmerList::receive_from_buffer(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count, int stage)  {
    
    auto& working_idx = working_idx_[src_process];
    auto end_idx = start_idx_[src_process+1];

    if (stage == 1) {
        end_idx = stage1_size_[src_process] + start_idx_[src_process];
    }

    size_t copy_num = std::min((max_count - current_count) / sizeof(KmerListEntryS) + 1, end_idx - working_idx);
    memcpy(kmerlist_.data() + working_idx, buffer + current_count, copy_num * sizeof(KmerListEntryS));
    current_count += copy_num * sizeof(KmerListEntryS);
    working_idx += copy_num;

    return working_idx == end_idx;
}

void GatheredKmerList::process(int nthr, int sort)  {
    size_t start_pos;
    sort_task(kmerlist_, sort, nthr, start_pos, std::accumulate(count_.begin(), count_.end(), 0));
    size_t valid_kmer;
    count_sorted_kmerlist(kmerlist_, filtered_kmerlist_, start_pos, std::accumulate(count_.begin(), count_.end(), 0), valid_kmer, true);

}





SendingGroup::SendingGroup(int destination, std::vector<std::shared_ptr<ScatteredTask>>& tasks) : 
    destination_(destination), tasks_(tasks)
{
    current_task_ = 0;

    ntasks_ = tasks.size();
}


bool SendingGroup::write_to_buffer(uint8_t* buffer, size_t max_count, int stage) {
    size_t current_count = 0;
    while (current_count <= max_count && current_task_ < ntasks_) {
        bool complete = tasks_[current_task_]->write_to_buffer(buffer, current_count, max_count, stage);
        if (complete) {
            current_task_++;
        }
    }

    return current_task_ == ntasks_;
}

void SendingGroup::init_exchange(int stage) {
    for (int i = 0; i < ntasks_; i++) {
        tasks_[i]->init_exchange(stage);
    }
    current_task_ = 0;
}



ReceivingGroup::ReceivingGroup(int source, std::vector<std::shared_ptr<GatheredTask>>& tasks) : 
    source_(source), tasks_(tasks)
{
    current_task_ = 0;
    // The task is responsible for keeping track of the sending info
}

void ReceivingGroup::receive_from_buffer(uint8_t* buffer, size_t max_count, int stage) {
    size_t current_count = 0;
    while (current_count <= max_count && current_task_ < tasks_.size()) {
        bool complete = tasks_[current_task_]->receive_from_buffer(buffer, source_, current_count, max_count, stage);
        if (complete) {
            current_task_++;
        }
    }

    return;
}

void ReceivingGroup::init_exchange(int stage) {
    if (source_ == 0){
        for (int i = 0; i < tasks_.size(); i++) {
            tasks_[i]->init_exchange(stage);
        }
    }
    current_task_ = 0;
}





TaskManager::TaskManager(TaskClassifier* classifier, TaskDispatcher* dispatcher, MPI_Comm comm, size_t batch_size, size_t max_send_size, int nthr) : 
        classifier_(classifier), dispatcher_(dispatcher), comm_(comm), logger_(comm),
        batch_size_(batch_size), max_send_size_(max_send_size), nthr_(nthr) {
    MPI_Comm_rank(comm, &myrank_);
    MPI_Comm_size(comm, &nprocs_);

    send_buffer1_ = new uint8_t[batch_size * nprocs_];
    send_buffer2_ = new uint8_t[batch_size * nprocs_];
    recv_buffer1_ = new uint8_t[batch_size * nprocs_];
    recv_buffer2_ = new uint8_t[batch_size * nprocs_];
}

TaskManager::~TaskManager() {
    delete classifier_;
    delete dispatcher_;
    delete[] send_buffer1_;
    delete[] send_buffer2_;
    delete[] recv_buffer1_;
    delete[] recv_buffer2_;
}

void TaskManager::init_scattered_tasks(std::vector<std::shared_ptr<ScatteredTask>>& sca_tasks) {
    sca_tasks_.clear();
    sca_tasks_ = sca_tasks;
    ntasks_ = sca_tasks.size();
}

void TaskManager::init_gathered_tasks(std::vector<std::shared_ptr<GatheredTask>>& gat_tasks) {
    gat_tasks_.clear();
    gat_tasks_ = gat_tasks;
}

void TaskManager::classify_tasks() {
    classifier_->classify(sca_tasks_, task_types_, comm_);

#if LOG_LEVEL >= 2
    Logger logger(comm_);
    for (int i = 0; i < task_types_.size(); i++) {
        logger() << "Task " << i << ": " << task_types_[i]<<std::endl;
    }
    logger.flush_root("Task classification result:");
#endif
}

void TaskManager::preprocess_tasks(int thr_per_worker) {

    size_t length_bytes = 0;
    size_t nsmtasks = 0;
    for (int i = 0; i < sca_tasks_.size(); i++) {
        if (task_types_[i] == 0) {
            nsmtasks++;
            auto supermer_task = std::dynamic_pointer_cast<ScatteredSupermers>(sca_tasks_[i]);
            length_bytes += supermer_task->get_size() * supermer_task->get_len_size();
        }
    }
    length_bytes /= nsmtasks;

    size_t stage1_round_est = length_bytes / max_send_size_;
    if (stage1_round_est == 0) {
        stage1_round_est = 1;
    }
    size_t stage1_kmerlist_cnt = stage1_round_est * max_send_size_ / sizeof(KmerListEntryS);
    MPI_Bcast(&stage1_kmerlist_cnt, 1, MPI_UNSIGNED_LONG_LONG, 0, comm_);

    stage1_kle_cnt_ = stage1_kmerlist_cnt;


    // TODO: parallelize this if there's performance issue
    for (int i = 0; i < sca_tasks_.size(); i++) {
        if (task_types_[i] == 0) {
            continue;
        } else if (task_types_[i] == 1) {
            auto supermer_task = std::dynamic_pointer_cast<ScatteredSupermers>(sca_tasks_[i]);
            sca_tasks_[i] = std::make_shared<ScatteredKmerList>(supermer_task, thr_per_worker, stage1_kmerlist_cnt);
            continue;
        }
        throw std::runtime_error("Unknown task type.");
    }
}

void TaskManager::dispatch() {
    dispatcher_->dispatch(comm_, sca_tasks_, task_destinations_);
    task_id_.clear();
    task_id_.resize(nprocs_);
    for (int i = 0; i < task_destinations_.size(); i++) {
        task_id_[task_destinations_[i]].push_back(i);
    }

#if LOG_LEVEL >= 2
    Logger logger(comm_);
    for (int i = 0; i < task_id_.size(); i++) {
        logger() << "Process " << i << ": ";
        for (int j = 0; j < task_id_[i].size(); j++) {
            logger() << task_id_[i][j] << " ";
        }
        logger()<<std::endl;
    }
    logger.flush_root("Task distribution:");
#endif
}

/// init the exchange groups, including the gathered tasks
void TaskManager::init_exchange_groups() {
    sending_groups_.clear();
    receiving_groups_.clear();

    auto& my_tasks = task_id_[myrank_];

    // exchange the count information
    std::vector<int32_t> scounts(nprocs_, 0);
    std::vector<int32_t> sdispls(nprocs_, 0);
    std::vector<int32_t> rcounts(nprocs_, 0);
    std::vector<int32_t> rdispls(nprocs_, 0);

    std::vector<size_t> sdata(ntasks_, 0);
    std::vector<size_t> rdata(nprocs_ * my_tasks.size(), 0);

    uint32_t offset = 0;
    for (int i = 0; i < nprocs_; i++) {
        scounts[i] = task_id_[i].size();
        sdispls[i] = offset;
        for (int j = 0; j < task_id_[i].size(); j++) {
            sdata[offset + j] = sca_tasks_[task_id_[i][j]]->get_size();
        }
        offset += task_id_[i].size();
    }

    for (int i = 0; i < nprocs_; i++) {
        rdispls[i] = i * my_tasks.size();
        rcounts[i] = my_tasks.size();
    }


    MPI_Alltoallv(sdata.data(), scounts.data(), sdispls.data(), MPI_UNSIGNED_LONG_LONG, 
            rdata.data(), rcounts.data(), rdispls.data(), MPI_UNSIGNED_LONG_LONG, comm_);

    std::vector<std::shared_ptr<GatheredTask>> gat_tasks;

    for (int i = 0; i < my_tasks.size(); i++) {
        std::vector<size_t> recv_counts;
        for (int j = 0; j < nprocs_; j++) {
            recv_counts.push_back(rdata[j * my_tasks.size() + i]);
        }

        auto taskid = task_id_[myrank_][i];
        if (task_types_[taskid] == 0) {
            gat_tasks.push_back(std::make_shared<GatheredSupermer>(taskid, 0, recv_counts));
        } else if (task_types_[taskid] == 1) {
            gat_tasks.push_back(std::make_shared<GatheredKmerList>(taskid, 1, recv_counts, stage1_kle_cnt_));
        }
    }

    init_gathered_tasks(gat_tasks);

    for (int i = 0; i < nprocs_; i++) {
        std::vector<std::shared_ptr<ScatteredTask>> tasks;
        for (int j = 0; j < task_id_[i].size(); j++) {
            tasks.push_back(sca_tasks_[task_id_[i][j]]);
        }
        sending_groups_.push_back(std::make_shared<SendingGroup>(i, tasks));
        receiving_groups_.push_back(std::make_shared<ReceivingGroup>(i, gat_tasks));
    }

}

void TaskManager::exchange(int stage) {
    exchange_init(stage);
    MPI_Barrier(comm_);
    int round = 0;
    while (true) {
        round++;

        if (exchange_progress(stage)) {
            break;
        }
    }
    exchange_clearup(stage);

    return;
}

void TaskManager::exchange_stage_finish() {
    delete[] classifier_;
    delete[] dispatcher_;
    delete[] send_buffer1_;
    delete[] send_buffer2_;
    delete[] recv_buffer1_;
    delete[] recv_buffer2_;

    classifier_ = nullptr;
    dispatcher_ = nullptr;
    send_buffer1_ = nullptr;
    send_buffer2_ = nullptr;
    recv_buffer1_ = nullptr;
    recv_buffer2_ = nullptr;

    sca_tasks_.clear();
    
    sending_groups_.clear();
    receiving_groups_.clear();
}

void TaskManager::process_tasks(int nworker, int nthr_tot, int sort) {
    int ntasks = gat_tasks_.size();

    int nworkers = std::min(nworker, ntasks);
    // int nworkers = nworker;
    int nthr = nthr_tot / nworkers;

    #pragma omp parallel for schedule(dynamic) num_threads(nworkers)
    for(int i = 0; i < ntasks; i++) {
        // TimerLocal timer;
        // timer.start();
        gat_tasks_[i]->process(nthr, sort);
        // logger_() << "Task " << i << " finished in time: " << timer.stop() << " on thread:" << omp_get_thread_num() <<std::endl;
    }

    // logger_.flush("Task processing:");

}

// get the size of gathered tasks' data in bytes
size_t TaskManager::get_gathered_task_bytes() {
    size_t sz = 0;
    for (int i = 0; i < gat_tasks_.size(); i++) {
        sz += gat_tasks_[i]->get_size_bytes();
    }
    return sz;
}

void TaskManager::copy_results(KmerListS& kmerlist) {
    size_t sz = 0;
    for (int i = 0; i < gat_tasks_.size(); i++) {
        sz += gat_tasks_[i]->get_kmerlist_size();
    }
    size_t offset = 0;
    kmerlist.resize(sz);

    for (int i = 0; i < gat_tasks_.size(); i++) {
        auto& task = gat_tasks_[i];
        auto& result = task->get_kmerlist();
        memcpy(kmerlist.data() + offset, result.data(), result.size() * sizeof(KmerListEntryS));
        offset += result.size();
    }

}

void TaskManager::exchange_init(int stage) {
    for(int i = 0; i < nprocs_; i++) {
        sending_groups_[i]->init_exchange(stage);
        receiving_groups_[i]->init_exchange(stage);
    }

    MPI_Barrier(comm_);
    
    write_sendbufs(stage);

    std::swap(send_buffer1_, send_buffer2_);
    std::swap(recv_buffer1_, recv_buffer2_);

    MPI_Ialltoall(send_buffer1_, batch_size_, MPI_BYTE, recv_buffer1_, batch_size_, MPI_BYTE, comm_, &req_);
}

bool TaskManager::exchange_progress(int stage) {

#if LOG_LEVEL >= 3
    Timer timer(comm_);
    timer.start();
#endif

    write_sendbufs(stage);

#if LOG_LEVEL >= 3
    timer.stop_and_log("Write sendbufs");
    timer.start();
#endif
#if LOG_LEVEL >= 4
    TimerLocal tm2;
    tm2.start();
#endif

    MPI_Wait(&req_, MPI_STATUS_IGNORE); //recved to 1 

#if LOG_LEVEL >= 4
    logger_() << "Round time: " << tm2.stop() << std::endl;
#endif
#if LOG_LEVEL >= 3
    timer.stop_and_log("MPI_Wait");
    timer.start();
#endif
#if LOG_LEVEL >= 4
    logger_.flush("Detailed time");
#endif

    std::swap(send_buffer1_, send_buffer2_);
    std::swap(recv_buffer1_, recv_buffer2_);

    MPI_Ialltoall(send_buffer1_, batch_size_, MPI_BYTE, recv_buffer1_, batch_size_, MPI_BYTE, comm_, &req_);

    bool ret = parse_recvbufs(stage);

#if LOG_LEVEL >= 3
    timer.stop_and_log("Parse recvbufs");
#endif

    return ret;
}

void TaskManager::exchange_clearup(int stage) {
    MPI_Wait(&req_, MPI_STATUS_IGNORE);
}

void TaskManager::write_sendbufs(int stage) {
    bool complete = true;

    // TODO: change back
    #pragma omp parallel for num_threads(nthr_)
    for(int i = 0; i < nprocs_; i++) {
        bool result = sending_groups_[i]->write_to_buffer(send_buffer2_ + batch_size_ * i, max_send_size_, stage);
        #pragma omp critical
        {
            complete = complete && result;
        }
    }

    for(int i = 0; i < nprocs_; i++) {
        if (complete) {
            send_buffer2_[batch_size_ * (i+1)-1] = 1;
        } else {
            send_buffer2_[batch_size_ * (i+1)-1] = 0;
        }
    }
}

bool TaskManager::parse_recvbufs(int stage) {
    bool complete = true;

    // TODO: change back
    #pragma omp parallel for num_threads(nthr_)
    for(int i = 0; i < nprocs_; i++) {
        receiving_groups_[i]->receive_from_buffer(recv_buffer2_ + batch_size_ * i, max_send_size_, stage);
        #pragma omp critical
        {
            complete = complete && *((bool*)(recv_buffer2_ + batch_size_ * (i+1)-1));
        }
    }

    return complete;
}


void FindKmerDestinationsParallel(const DnaBuffer& myreads, int nthreads, int tot_tasks, ParallelData& data) {
    assert(nthreads > 0);

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (size_t i = 0; i < myreads.size(); ++i) {
        int tid = omp_get_thread_num();

        /* register the read for storage */
        auto &dest = data.register_new_read(tid, i);
        if (myreads[i].size() < KMER_SIZE)
            continue;
        dest.reserve(myreads[i].size() - KMER_SIZE + 1);

        /* prepare the vars */
        std::vector<TMmer> repmers = TMmer::GetRepMmers(myreads[i]);
        Minimizer_Deque deque;
        int head_pos = 0;

        /* initialize the deque */
        for(; head_pos < KMER_SIZE - MINIMIZER_SIZE; head_pos++) {
            deque.insert_minimizer(repmers[head_pos].GetHash(), head_pos);
        }
        int tail_pos = head_pos - KMER_SIZE + MINIMIZER_SIZE - 1;

        /* start recording the destinations */
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

void Minimizer_Deque::remove_minimizer(int pos) {
    while (!deq.empty() && deq.front().second <= pos) {
        deq.pop_front();
    }
}

void Minimizer_Deque::insert_minimizer(uint64_t hash, int pos) {
    while (!deq.empty() && deq.back().first > hash) {
        deq.pop_back();
    }
    deq.push_back(std::make_pair(hash, pos));
}

uint64_t Minimizer_Deque::get_current_minimizer() {
    return deq.front().first;
}


ParallelData::ParallelData(int nthr) {
    this->nthr = nthr;
    destinations.resize(nthr);
    readids.resize(nthr);
}

std::vector<std::vector<int>>& ParallelData::get_my_destinations(int tid) {
    return destinations[tid];
}

std::vector<uint64_t>& ParallelData::get_my_readids(int tid) {
    return readids[tid];
}

std::vector<int>& ParallelData::register_new_read(int tid, uint64_t readid) {
    destinations[tid].push_back(std::vector<int>());
    readids[tid].push_back(readid);
    return destinations[tid].back();
}

void SupermerEncoder::copy_bits(std::vector<uint8_t>& dst, const uint8_t* src, uint64_t start_pos, int len){
    size_t start = dst.size();
    dst.resize(start + cnt_bytes(len)); // TODO: reserve beforehand if there's a performance issue

    for (int i = 0; i < len; i++) {
        /* the order is confusing, just make sure we keep it consistent */
        int loc = i + start_pos;
        bool first_bit = (src[loc / 4] >> (7 - 2*(loc%4))) & 1;
        bool second_bit = (src[loc / 4] >> (6 - 2*(loc%4))) & 1;
        dst[start + i/4] |= (first_bit << (7 - (i%4)*2)) | (second_bit << (6 - (i%4)*2));
    }
}

void SupermerEncoder::encode(const std::vector<int>& dest, const DnaSeq& read){
    if (read.size() < KMER_SIZE) return;

    /* initial conditions */
    uint32_t start_pos = 0;
    int cnt = 1;
    int last_dst = dest[0];

    /* iterate through the destinations */
    for (int i = 1; i <= dest.size(); i++) {

        if(i == dest.size() || dest[i] != last_dst || cnt == max_supermer_len_ - KMER_SIZE + 1) {
            /* encode the supermer */
            size_t len = cnt + KMER_SIZE - 1;
            (*tasks_[last_dst])
                    .get_length_buffer(tid_)
                    .push_back(len);
            auto& supermer_buffer = (*tasks_[last_dst])
                    .get_supermer_buffer(tid_);
            copy_bits(supermer_buffer, read.data(), start_pos, len);

            /* reset the counter */
            if (i < dest.size()) last_dst = dest[i];
            cnt = 0;
            start_pos = i;
        }

        /* increment the counter */
        cnt++;
    }
}

void PlainClassifier::classify(std::vector<std::shared_ptr<ScatteredTask>> sca_tasks, std::vector<int>& task_types, MPI_Comm comm)  {
    task_types.resize(sca_tasks.size(), 0);
    for (int i = 0; i < sca_tasks.size(); i++) {
        task_types[i] = 0;
    }
}

void HeavyHitterClassifier::classify(std::vector<std::shared_ptr<ScatteredTask>> sca_tasks, std::vector<int>& task_types, MPI_Comm comm)  {
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int ntask = sca_tasks.size();
    task_types.resize(ntask, 0);

    std::vector<size_t> task_sizes(ntask, 0);
    for (int i = 0; i < sca_tasks.size(); i++) {
        task_sizes[i] = std::reinterpret_pointer_cast<ScatteredSupermers>(sca_tasks[i])->get_size_kmer();
    }

    std::vector<size_t> task_sizes_global(ntask, 0);
    MPI_Reduce(task_sizes.data(), task_sizes_global.data(), ntask, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);

    if (myrank == 0) {
        task_sizes = task_sizes_global;
        size_t total_size = 0;
        for(int i = 0; i < ntask; i++){
            total_size += task_sizes[i];
        }
        size_t avg_size = total_size / ntask;
        for (int i = 0; i < ntask; i++) {
            if (task_sizes[i] > avg_size * UNBALANCED_RATIO) {
                task_types[i] = 1;
            } else {
                task_types[i] = 0;
            }
        }
    }

    MPI_Bcast(task_types.data(), ntask, MPI_INT, 0, comm);
}


bool BalancedDispatcher::try_dispatch(std::vector<TaskInfo>& task_info, std::vector<int>& destinations, int nprocs, size_t avg_size, double coe) {
    int ntasks = task_info.size();
    size_t upper_bound = avg_size * coe;
    destinations.resize(ntasks, -1);


    std::vector<std::vector<int>> assigned_task_id(nprocs);
    std::vector<size_t> assigned_size(nprocs, 0);
    std::vector<bool> assigned(ntasks, false);

    auto ppn = get_ppn();
    auto nodes = nprocs / ppn;
    /* stage1: each process get assigned a big task */
    for(int i = 0; i < nprocs; i++) {
        //auto rank = (i % nodes) * ppn + (i / nodes + i) % ppn;
        auto rank = i;
        int id = ntasks - 1 - i;
        assigned[id] = true;
        assigned_task_id[rank].push_back(task_info[id].id); // change here (ver a)
        assigned_size[rank] += task_info[id].sz * task_info[id].coe;
    }

    /* stage2: assigned the result of the tasks */
    int cur_proc = nprocs - 1;
    for(int i = 0; i < ntasks; i++) {
        if (assigned[i]) { continue; }

        int cnt = 0;
        while (cnt < nprocs) {
            if (assigned_size[cur_proc] + task_info[i].sz <= upper_bound) {
                assigned_task_id[cur_proc].push_back(task_info[i].id);
                assigned_size[cur_proc] += task_info[i].sz * task_info[i].coe;
                assigned[i] = true;
                cur_proc--;
                if (cur_proc < 0) {
                    cur_proc += nprocs;
                }
                break;
            } else {
                cur_proc--;
                if (cur_proc < 0) {
                    cur_proc += nprocs;
                }
                cnt++;
            }
        }
        if (cnt == nprocs) {
            return false;
        }
    }

    /* stage3: we're sure all tasks are assigned if reach here */
    for (int i = 0; i < nprocs; i++) {
        for (int j = 0; j < assigned_task_id[i].size(); j++) {
            destinations[assigned_task_id[i][j]] = i;
        }
    }
    return true;
}

void BalancedDispatcher::dispatch(MPI_Comm comm, std::vector<std::shared_ptr<ScatteredTask>>& tasks, std::vector<int>& destinations)  {
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
    
    /* gather the task size information to rank 0 */
    std::vector<size_t> task_size;
    for (int i = 0; i < tasks.size(); i++) {
        task_size.push_back(tasks[i]->get_size_bytes());
    }

    std::vector<size_t> task_size_global(task_size.size(), 0);
    MPI_Reduce(task_size.data(), task_size_global.data(), task_size.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    
    /* rank 0 dispatch the tasks */
    if (myrank == 0) {
        task_size = task_size_global;
        for(int i = 0; i < task_size.size(); i++) {
            std::cout<<"Task "<<i<<" size: "<<task_size[i]<<std::endl;
        }

        std::vector<TaskInfo> task_info;
        for (int i = 0; i < tasks.size(); i++) {
            task_info.push_back({i, task_size[i], tasks[i]->get_extra_coe()});
        }
        std::sort(task_info.begin(), task_info.end());

        size_t total_size = 0;
        for (int i = 0; i < task_size.size(); i++) {
            total_size += task_size[i];
        }
        size_t avg_size = total_size / nprocs;

        double coe = 1.0;
        bool success = false;
        while( coe < DISPATCH_UPPER_COE ) {
            coe += DISPATCH_STEP;
            if (try_dispatch(task_info, destinations, nprocs, avg_size, coe)) {
                success = true;
                break;
            }
        }
        if (!success) {
            throw std::runtime_error("Cannot dispatch tasks. May be too unbalanced.");
        }
    }

    destinations.resize(tasks.size(), -1);

    MPI_Bcast(destinations.data(), tasks.size(), MPI_INT, 0, comm);

}

void print_kmer_histogram(const KmerListS& kmerlist, MPI_Comm comm) {
    #if LOG_LEVEL >= 2

    Logger logger(comm);
    size_t maxcount = std::accumulate(kmerlist.cbegin(), kmerlist.cend(), 0, [](size_t cur, const auto& entry) { return std::max(cur, entry.cnt); });

    MPI_Allreduce(MPI_IN_PLACE, &maxcount, 1, MPI_INT, MPI_MAX, comm);

    std::vector<int> histo(maxcount+1, 0);

    for(size_t i = 0; i < kmerlist.size(); ++i)
    {
        int cnt = kmerlist[i].cnt;
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


int sort_decision(size_t total_bytes, Logger& logger) {
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
            if (memfree >= total_bytes) {
                sort = 2;
                logger() << "Enough memory available. Using RADULS for sorting." ;
            } else {
                logger() << "Not enough memory available. Using PARADIS for sorting.";
                sort = 1;
            }
        }
    }
    return sort;
}


template<typename T>
void sort_task(std::vector<T>& kmerseeds, int sort, int thr_per_worker, size_t& start_pos, size_t seedcnt) {
    if (sort == 1){
        start_pos= 0;
        paradis::sort<T, TKmer::NBYTES>(kmerseeds.data(), kmerseeds.data() + seedcnt, thr_per_worker);
    } else {
        uint8_t* tmp_arr = new uint8_t[seedcnt * sizeof(T) + 256];
        uint8_t* tmp = tmp_arr + 256 - (size_t)tmp_arr % 256;
        raduls::CleanTmpArray(tmp, seedcnt, sizeof(T), thr_per_worker);

        // std::cout<<"z1"<<std::endl;
        // RADULS needs padding. reserved in bucket assistant
        uint8_t* start = (uint8_t*)kmerseeds.data();
        // std::cout<<"start is"<<(size_t)start<<std::endl;
        int cnt = 0;
        while( (size_t)start % 256 != 0) {
            start += sizeof(T);
            kmerseeds.push_back(kmerseeds[cnt]);
            cnt++;
        }
        start_pos= cnt;

        // std::cout<<"z2 "<<sizeof(T)<<" "<<TKmer::NBYTES<<" "<<seedcnt<<" "<<thr_per_worker<<std::endl;

        // std::cout<<"some data"<<((size_t)start)<< " " << ((size_t)tmp)<< " " << kmerseeds[start_pos].kmer << std::endl;

        // std::cout<<kmerseeds.size()<<" "<<seedcnt<<std::endl;

        raduls::RadixSortMSD(start, tmp, seedcnt, sizeof(T), TKmer::NBYTES, thr_per_worker);

        // std::cout<<"pass here"<<((size_t)start)<< " " << ((size_t)tmp)<< " " << kmerseeds[start_pos].kmer << std::endl;

        
        delete[] (tmp_arr);
    }

}


void count_sorted_kmers(std::vector<KmerSeedStruct>& kmerseeds, KmerListS& kmerlist, size_t start_pos, size_t seedcnt, size_t& valid_kmer, bool filter) {
    kmerlist.clear();
    kmerlist.reserve(seedcnt / LOWER_KMER_FREQ);
    valid_kmer = 0;

    TKmer last_mer = kmerseeds[start_pos].kmer;
    uint64_t cur_kmer_cnt = 1;
    for (size_t idx = start_pos + 1; idx <= start_pos + seedcnt; idx++) {
        TKmer cur_mer;
        if (idx != start_pos + seedcnt) {
            cur_mer = kmerseeds[idx].kmer;
        }
        
        if (cur_mer == last_mer && idx != start_pos + seedcnt) {
            cur_kmer_cnt++;
            continue;
        } 

        //std::cerr<<"cur_kmer_cnt "<<cur_kmer_cnt<<std::endl;
        if (!filter || (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) ) {
            kmerlist.emplace_back(last_mer, cur_kmer_cnt);
            valid_kmer++;
        }

        cur_kmer_cnt = 1;
        last_mer = cur_mer;
    }
}

void count_sorted_kmerlist(KmerListS& kmers, KmerListS& kmerlist, size_t start_pos, size_t seedcnt, size_t& valid_kmer, bool filter) {
    kmerlist.clear();
    kmerlist.reserve(seedcnt / LOWER_KMER_FREQ);
    valid_kmer = 0;

    TKmer last_mer = kmers[start_pos].kmer;
    uint64_t cur_kmer_cnt = kmers[start_pos].cnt;
    for (size_t idx = start_pos + 1; idx <= start_pos + seedcnt; idx++) {
        TKmer cur_mer;
        if (idx != start_pos + seedcnt) {
            cur_mer = kmers[idx].kmer;
        }
        
        if (cur_mer == last_mer && idx != start_pos + seedcnt) {
            cur_kmer_cnt += kmers[idx].cnt;
            continue;
        } 

        // if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ)
        if (!filter || (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) ) {
            kmerlist.emplace_back(last_mer, cur_kmer_cnt);
            valid_kmer++;
        }

        cur_kmer_cnt = kmers[idx].cnt;
        last_mer = cur_mer;
    }
}