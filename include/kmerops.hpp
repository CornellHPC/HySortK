#ifndef HYSORTK_KMEROPS_H_
#define HYSORTK_KMEROPS_H_

#include <mpi.h>
#include <omp.h>
#include <deque>
#include <cmath>
#include "kmer.hpp"
#include "timer.hpp"
#include "dnaseq.hpp"
#include "dnabuffer.hpp"
#include "logger.hpp"
#include "compiletime.h"

#include "supermer.hpp"

namespace hysortk {

/*
 * Constants
 */
#define PAD_SIZE 256        // For raduls

/* 
 * Assistant Functions
 */

/// @brief Assistance function to get the number of processes per node
int get_ppn();


/// @brief Assistance function to calculate how many dna base to pad
inline int pad_base(const int& len) {
    return (4 - len % 4);
}


/// @brief Assistance function to calculate how many bytes
inline int cnt_bytes(const int& len) {
    return (len + pad_base(len)) / 4;
}


/*
 * Supermer Functions
 */


/// @brief Assistant function to assign a minimizer a task
int GetMinimizerOwner(const uint64_t& hash, int tot_tasks);


/// @brief Not in use. Assistant function to assign a kmer a task
int GetKmerOwner(const TKmer& kmer, int nprocs);


/// @brief Assistant class to store the data from different threads when parsing the reads
class ParallelData{
private:
    int nthr;
    std::vector<std::vector<std::vector<int>>> destinations;
    std::vector<std::vector<ReadId>> readids;
public:
    ParallelData(int nthr);

    std::vector<std::vector<int>>& get_my_destinations(int tid);
    std::vector<ReadId>& get_my_readids(int tid);
    std::vector<int>& register_new_read (int tid, ReadId readid);

};


/// @brief Function that parse the reads 
void FindKmerDestinationsParallel(const DnaBuffer& myreads, int nthreads, int tot_tasks, ParallelData& data);


/// @brief Assistant class to help determine the minimizers of kmers from a read
struct Minimizer_Deque {
    /* The first item is the hash value, the second one is the position of the minimizer */
    std::deque<std::pair<uint64_t, int>> deq;

    void remove_minimizer(int pos);
    void insert_minimizer(uint64_t hash, int pos);
    uint64_t get_current_minimizer();
};


/*
 * Task-related classes and functions
 */

class ScatteredTask {
protected:
    int taskid_;
    int tasktype_;
    int destination_;

    ScatteredTask(int taskid, int tasktype) : taskid_(taskid), tasktype_(tasktype) {}
public:
    virtual size_t get_size_bytes() = 0;
    virtual size_t get_size() = 0;
    virtual int get_extra_coe() {return 1;};
    virtual bool write_to_buffer(uint8_t* buffer, size_t& current_count, size_t max_count, int stage) = 0;
    virtual void init_exchange(int stage) = 0;
    int get_type() {
        return tasktype_;
    }

    int get_taskid() {
        return taskid_;
    }
};

class ScatteredSupermers : public ScatteredTask {
private:
    int nthr_;
    int extra_coe_ = 1;

    std::vector<std::vector<uint8_t>> supermers_;
    std::vector<std::vector<length_t>> lengths_;
    constexpr static int LEN_SIZE = sizeof(length_t);

    int working_thr_;
    int working_idx_;
    int working_buf_idx_;

    bool write_to_buffer_stage1(uint8_t* buffer, size_t &current_count, size_t max_count);
    bool write_to_buffer_stage2(uint8_t* buffer, size_t &current_count, size_t max_count);

public:
    ScatteredSupermers(int taskid, int nthr);

    int get_nthr() { return nthr_; }
    std::vector<uint8_t>& get_supermer_buffer (int tid) { return supermers_[tid]; }
    std::vector<length_t>& get_length_buffer (int tid) { return lengths_[tid]; }   

    // This function is getting the number of supermers, not the number of kmers
    size_t get_size_bytes() override;
    size_t get_size() override;
    size_t get_len_size() { return LEN_SIZE; }
    size_t get_size_kmer();
    int get_extra_coe() override { return extra_coe_; }
    void init_exchange(int stage);
    bool write_to_buffer(uint8_t* buffer, size_t& current_count, size_t max_count, int stage) override;
};



class ScatteredKmerList : public ScatteredTask {
private:
    KmerListS kmerlist_;
    size_t working_idx_;

    size_t stage1_size_;

public:
    ScatteredKmerList(std::shared_ptr<ScatteredSupermers> supermer_task, int thr_per_worker=4, size_t stage1_size=0);
    
    size_t get_size_bytes() override { return kmerlist_.size() * sizeof(KmerListEntryS); }
    size_t get_size() override { return kmerlist_.size(); }
    void init_exchange(int stage) override { return; }
    int get_extra_coe() override { return 1; }
    bool write_to_buffer(uint8_t* buffer, size_t& current_count, size_t max_count, int stage) override;
};


class GatheredTask {
protected:
    int taskid_;
    int tasktype_;
    std::vector<size_t> count_;

    std::vector<size_t> start_idx_;
    std::vector<size_t> working_idx_;

    KmerListS filtered_kmerlist_;


public:
    GatheredTask(int taskid, int tasktype, std::vector<size_t> count);

    virtual void init_exchange(int stage) = 0;
    virtual bool receive_from_buffer(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count, int stage) = 0;
    virtual void process(int nthr, int sort) = 0;
    virtual size_t get_size_bytes() = 0;
    size_t get_kmerlist_size() { return filtered_kmerlist_.size(); }
    KmerListS& get_kmerlist() { return filtered_kmerlist_;}
};


class GatheredSupermer : public GatheredTask {
private:
    std::vector<KmerSeedStruct> kmer_;
    std::vector<length_t> length_;
    std::vector<size_t> kmer_count_;
    std::vector<size_t> working_kmer_idx_;

    constexpr static int LEN_SIZE = sizeof(length_t);

    bool receive_from_buffer_stage1(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count);
    bool receive_from_buffer_stage2(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count);

public:
    GatheredSupermer(int taskid, int tasktype, std::vector<size_t> count) : GatheredTask(taskid, tasktype, count) {}
    void init_exchange(int stage) override;
    bool receive_from_buffer(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count, int stage) override;
    void process(int nthr, int sort) override ;
    size_t get_size_bytes() override { return std::accumulate(kmer_count_.begin(), kmer_count_.end(), 0) * sizeof(KmerSeedStruct); }

};



class GatheredKmerList : public GatheredTask {
private:
    KmerListS kmerlist_;
    std::vector<size_t> stage1_size_;

public:
    GatheredKmerList(int taskid, int tasktype, std::vector<size_t> count, size_t stage1_size=0);
    void init_exchange(int stage) override { return; }
    bool receive_from_buffer(uint8_t* buffer, int src_process, size_t &current_count, size_t max_count, int stage) override;
    void process(int nthr, int sort) override;
    size_t get_size_bytes() override { return std::accumulate(count_.begin(), count_.end(), 0) * sizeof(KmerListEntryS); }
};


/// @brief Assistant class to encode the supermer to bytes
struct SupermerEncoder{
    std::vector<std::shared_ptr<ScatteredSupermers>> tasks_;
    int max_supermer_len_;
    int tid_;

    SupermerEncoder(std::vector<std::shared_ptr<ScatteredTask>> tasks, 
            int tid, 
            int max_supermer_len) : 
        tid_(tid), max_supermer_len_(max_supermer_len) {
            for (auto& task : tasks) {
                tasks_.push_back(std::dynamic_pointer_cast<ScatteredSupermers>(task));
            }
        }

    /* no practical solution for storing bits instead of bytes is found */
    void copy_bits(std::vector<uint8_t>& dst, const uint8_t* src, uint64_t start_pos, int len);

    /* encode the supermer to the buffer */
    void encode(const std::vector<int>& dest, const DnaSeq& read, ReadId readid);
};


/// @brief Base class to classify the tasks
class TaskClassifier{
public:
    virtual void classify(std::vector<std::shared_ptr<ScatteredTask>> sca_tasks, std::vector<int>& task_types, MPI_Comm comm) = 0;
};

/// @brief All tasks are considered the plain task
class PlainClassifier : public TaskClassifier {
public:
    void classify(std::vector<std::shared_ptr<ScatteredTask>> sca_tasks, std::vector<int>& task_types, MPI_Comm comm) override;
};

/// @brief Classify the tasks based on the size
class HeavyHitterClassifier : public TaskClassifier {
public:
    void classify(std::vector<std::shared_ptr<ScatteredTask>> sca_tasks, std::vector<int>& task_types, MPI_Comm comm) override;
};

/// @brief Base class that provides a strategy to dispatch the tasks to MPI processes
class TaskDispatcher {
public:
    virtual void dispatch(MPI_Comm comm, std::vector<std::shared_ptr<ScatteredTask>>& tasks, std::vector<int>& destinations) = 0;
};

/// @brief Dispatcher that will dispatch the tasks in a round-robin way
class RoundRobinDispatcher : public TaskDispatcher {
public:
    void dispatch(MPI_Comm comm, std::vector<std::shared_ptr<ScatteredTask>>& tasks, std::vector<int>& destinations) override;
};

/// @brief Dispatcher that will dispatch the tasks in a balanced way
class BalancedDispatcher : public TaskDispatcher {
private:
    struct TaskInfo{
        int64_t id;
        size_t sz;
        int64_t coe;
        bool operator < (const TaskInfo& o) const { return sz < o.sz; }
        bool operator > (const TaskInfo& o) const { return sz > o.sz; }
        bool operator == (const TaskInfo& o) const { return sz == o.sz; }
    };

    /// @brief Try to dispatch the tasks with a given upper bound
    bool try_dispatch(std::vector<TaskInfo>& task_info, std::vector<int>& destinations, int nprocs, size_t avg_size, double coe);

public:
    void dispatch(MPI_Comm comm, std::vector<std::shared_ptr<ScatteredTask>>& tasks, std::vector<int>& destinations) override;
};


class SendingGroup {
private:
    int destination_;
    std::vector<std::shared_ptr<ScatteredTask>> tasks_;
    int ntasks_;
    int current_task_;

public:
    SendingGroup(int destination, std::vector<std::shared_ptr<ScatteredTask>>& tasks);
    bool write_to_buffer(uint8_t* buffer, size_t max_count, int stage);
    void init_exchange(int stage);
};



class ReceivingGroup {
private:
    int source_;
    std::vector<std::shared_ptr<GatheredTask>> tasks_;
    int current_task_;

public:
    ReceivingGroup(int source, std::vector<std::shared_ptr<GatheredTask>>& tasks);
    void receive_from_buffer(uint8_t* buffer, size_t max_count, int stage);
    void init_exchange(int stage);
};



class TaskManager {
private:
    MPI_Comm comm_;
    Logger logger_;
    int myrank_;
    int nprocs_;
    int ntasks_;

    TaskClassifier* classifier_;
    TaskDispatcher* dispatcher_;

    size_t batch_size_;
    size_t max_send_size_;
    size_t stage1_kle_cnt_;            // How many kmerlist entry is allowed to send in stage 1
    int nthr_;
    uint8_t* send_buffer1_ = nullptr;
    uint8_t* send_buffer2_ = nullptr;
    uint8_t* recv_buffer1_ = nullptr;
    uint8_t* recv_buffer2_ = nullptr;
    MPI_Request req_;


    std::vector<std::shared_ptr<ScatteredTask>> sca_tasks_;
    std::vector<std::shared_ptr<GatheredTask>> gat_tasks_;
    std::vector<int> task_types_;
    std::vector<int> task_destinations_;
    std::vector<std::vector<int>> task_id_;

    std::vector<std::shared_ptr<SendingGroup>> sending_groups_;
    std::vector<std::shared_ptr<ReceivingGroup>> receiving_groups_;



public:
    TaskManager(TaskClassifier* classifier, TaskDispatcher* dispatcher, MPI_Comm comm, size_t batch_size, size_t max_send_size, int nthr);
    ~TaskManager();

    void init_scattered_tasks(std::vector<std::shared_ptr<ScatteredTask>>& sca_tasks);
    void init_gathered_tasks(std::vector<std::shared_ptr<GatheredTask>>& gat_tasks);
    void classify_tasks();
    void preprocess_tasks(int thr_per_worker=4);
    void dispatch();
    /// init the exchange groups, including the gathered tasks
    void init_exchange_groups();
    void exchange(int stage);
    void process_tasks(int nworker, int nthr_tot, int sort) ;

    // get the size of gathered tasks' data in bytes
    size_t get_gathered_task_bytes();
    void copy_results(KmerListS& kmerlist);

    void exchange_stage_finish();

private:
    // exchange related functions
    void exchange_init(int stage);
    bool exchange_progress(int stage);
    void exchange_clearup(int stage);
    void write_sendbufs(int stage);
    bool parse_recvbufs(int stage);
};


/// @brief Decide the sort algorithm based on the memory available
int sort_decision(size_t total_bytes, Logger& logger);

/// @brief Sort the kmer seeds
template<typename T>
void sort_task(std::vector<T>& kmerseeds, int sort, int thr_per_worker, size_t& start_pos, size_t seedcnt);

/// @brief Count the sorted kmer seeds
void count_sorted_kmers(std::vector<KmerSeedStruct>& kmerseeds, KmerListS& kmerlist, size_t start_pos, size_t seedcnt, size_t& valid_kmer, bool filter=true);

/// @brief Count the sorted kmer lists
void count_sorted_kmerlist(KmerListS& kmers, KmerListS& kmerlist, size_t start_pos, size_t seedcnt, size_t& valid_kmer, bool filter=true);

std::shared_ptr<TaskManager> prepare_supermer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int avg_thr_per_worker=THREAD_PER_WORKER,
     int max_thr_m=MAX_THREAD_MEMORY_BOUNDED);

void exchange_supermer(std::shared_ptr<TaskManager> task_manager, MPI_Comm comm, int avg_thr_per_worker=THREAD_PER_WORKER);

std::unique_ptr<KmerListS> filter_kmer(std::shared_ptr<TaskManager> task_manager, MPI_Comm comm, int avg_thr_per_worker=THREAD_PER_WORKER);

} // namespace hysortk

#endif // KMEROPS_H_