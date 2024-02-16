#ifndef KMEROPS_H_
#define KMEROPS_H_

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

typedef uint32_t PosInRead;
typedef  int64_t ReadId;

typedef std::array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef std::array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef std::tuple<TKmer, int> KmerListEntry;
typedef std::vector<KmerListEntry> KmerList;

struct KmerSeedStruct{
    TKmer kmer;      

    KmerSeedStruct(TKmer kmer) : kmer(kmer) {};
    KmerSeedStruct(const KmerSeedStruct& o) : kmer(o.kmer) {};
    KmerSeedStruct(KmerSeedStruct&& o) :  
        kmer(std::move(o.kmer)) {};
    KmerSeedStruct() {};

    int GetByte(int &i) const { return kmer.getByte(i); }
    bool operator < (const KmerSeedStruct& o) const { return kmer < o.kmer; }
    bool operator == (const KmerSeedStruct& o) const { return kmer == o.kmer; }
    bool operator != (const KmerSeedStruct& o) const { return kmer != o.kmer; }

    KmerSeedStruct& operator=(const KmerSeedStruct& o)
    {
        kmer = o.kmer;
        return *this;
    }
};

typedef std::vector<std::vector<KmerSeedStruct>> KmerSeedBuckets;
typedef std::vector<std::vector<std::vector<KmerSeedStruct>>> KmerSeedVecs;

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = THREAD_PER_TASK);

int GetKmerOwner(const TKmer& kmer, int nprocs);

void print_kmer_histogram(const KmerList& kmerlist, MPI_Comm comm);

class ParallelData{
    private:

        std::vector<std::vector<std::vector<int>>> destinations;
        std::vector<std::vector<uint64_t>> readids;

    public:
        int nprocs;
        int ntasks;
        int nthr_membounded;
        std::vector<std::vector<std::vector<uint32_t>>> lengths;
        std::vector<std::vector<std::vector<uint8_t>>> supermers;

        ParallelData(int nprocs, int ntasks, int nthr_membounded) {
            this->nprocs = nprocs;
            this->ntasks = ntasks;
            this->nthr_membounded = nthr_membounded;
            destinations.resize(nthr_membounded);
            lengths.resize(nthr_membounded);
            supermers.resize(nthr_membounded);
            readids.resize(nthr_membounded);

            for (int i = 0; i < nthr_membounded; i++) {
                lengths[i].resize(nprocs * ntasks);
                supermers[i].resize(nprocs * ntasks);
            }
        }

        std::vector<std::vector<int>>& get_my_destinations(int tid) {
            return destinations[tid];
        }

        std::vector<int>& register_new_destination (int tid, uint64_t readid) {
            destinations[tid].push_back(std::vector<int>());
            readids[tid].push_back(readid);
            return destinations[tid].back();
        }

        std::vector<std::vector<uint32_t>>& get_my_lengths(int tid) {
            return lengths[tid];
        }

        std::vector<std::vector<uint8_t>>& get_my_supermers(int tid) {
            return supermers[tid];
        }

        std::vector<uint64_t>& get_my_readids(int tid) {
            return readids[tid];
        }

        size_t get_supermer_cnt(int procid, int taskid) {
            size_t cnt = 0;
            for (int i = 0; i < nthr_membounded; i++) {
                cnt += lengths[i][procid * ntasks + taskid].size();
            }
            return cnt;
        }
};

ParallelData
prepare_supermer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int thr_per_task = THREAD_PER_TASK,
     int max_thr_membounded = MAX_THREAD_MEMORY_BOUNDED);

inline int pad_bytes(const int& len) {
    return (4 - len % 4) * 2;
}

inline int cnt_bytes(const int& len) {
    return (len + pad_bytes(len)) / 4;
}

struct SupermerEncoder{
    std::vector<std::vector<uint32_t>>& lengths;
    std::vector<std::vector<uint8_t>>& supermers;
    int max_supermer_len;

    SupermerEncoder(std::vector<std::vector<uint32_t>>& lengths, 
            std::vector<std::vector<uint8_t>>& supermers, 
            int max_supermer_len) : 
            lengths(lengths), supermers(supermers), max_supermer_len(max_supermer_len) {};


    // std::vector<bool> is depreciated, std::deque<bool> does not guarantee contiguity
    // std::bitset has a fixed length, so we use std::vector<uint8_t> instead
    // maybe we should skip this part and use a bitset for sendbuf directly
    void copy_bits(std::vector<uint8_t>& dst, const uint8_t* src, uint64_t start_pos, int len){

        size_t start = dst.size();
        dst.resize(start + cnt_bytes(len));

        for (int i = 0; i < len; i++) {
            /* the order is confusing, just make sure we keep it consistent*/
            int loc = i + start_pos;
            bool first_bit = (src[loc / 4] >> (7 - 2*(loc%4))) & 1;
            bool second_bit = (src[loc / 4] >> (6 - 2*(loc%4))) & 1;
            dst[start + i/4] |= (first_bit << (7 - (i%4)*2)) | (second_bit << (6 - (i%4)*2));
        }

    }


    void encode(const std::vector<int>& dest, const DnaSeq& read){
        
        /* initial conditions */
        uint32_t start_pos = 0;
        int cnt = 1;
        int last_dst = dest[0];

        for (int i = 1; i <= dest.size(); i++) {

            if(i == dest.size() || dest[i] != last_dst || cnt == max_supermer_len - KMER_SIZE + 1) {
                /* encode the supermer */
                size_t len = cnt + KMER_SIZE - 1;
                lengths[last_dst].push_back(len);
                copy_bits(supermers[last_dst], read.data(), start_pos, len);

                /* reset the counter */
                
                if (i < dest.size()) last_dst = dest[i];
                cnt = 0;
                start_pos = i;
            }

            /* increment the counter */
            cnt++;
        }
    }
};

/* this is designed to be a universal all to all batch exchanger which supports group */
class BatchExchanger
{
protected:
    MPI_Comm comm;
    MPI_Request req;

    int nprocs;
    int myrank;
    int ntasks;                     /* number of tasks */
    int round;
    double t1, t2, t3, t4;

    size_t batch_size;              /* maxium batch size in bytes */
    size_t send_limit;              /* maxium size of meaningful data in bytes */
    size_t max_element_size;        /* maxium size of a single element in bytes */

    uint8_t* sendbuf_x;
    uint8_t* recvbuf_x;
    uint8_t* sendbuf_y;
    uint8_t* recvbuf_y;

    virtual bool write_sendbuf(uint8_t* addr, int procid) = 0;
    virtual void parse_recvbuf(uint8_t* addr, int procid) = 0;

    void write_sendbufs(uint8_t* addr) {
        bool complete = true;
        for (int i = 0; i < nprocs; i++) {
            bool this_complete = write_sendbuf(addr + i * batch_size, i);
            if (!this_complete) {
                complete = false;
            }
        }

        for(int i = 0; i < nprocs; i++) {
            addr[(i+1) * batch_size - 1] = complete ? 1 : 0;
        } 
    }

    virtual void parse_recvbufs(uint8_t* addr) {
        for (int i = 0; i < nprocs; i++) {
            parse_recvbuf(addr + i * batch_size, i);
        }
    }

    bool check_complete(uint8_t* recvbuf) {
        bool flag = true;
        for (int i = 0; i < nprocs; i++) {
            if (recvbuf[(i+1) * batch_size - 1] == 0) {
                flag = false;
                break;
            }
        }
        return flag;
    }

public:
    enum Status
    {
        BATCH_NOT_INIT = 0,
        BATCH_SENDING = 1,
        BATCH_DONE = 2
    } status;


    void initialize() {
        if (status != BATCH_NOT_INIT) {
            return;
        }

        sendbuf_x = new uint8_t[batch_size * nprocs];
        recvbuf_x = new uint8_t[batch_size * nprocs];
        sendbuf_y = new uint8_t[batch_size * nprocs];
        recvbuf_y = new uint8_t[batch_size * nprocs];

        status = BATCH_SENDING;
        write_sendbufs(sendbuf_y);
        MPI_Ialltoall(sendbuf_y, batch_size, MPI_BYTE, recvbuf_y, batch_size, MPI_BYTE, comm, &req);
    }


    void progress() {
        if (status != BATCH_SENDING) {
            return;
        }
        round++;

        /* x is the buffer we can use to write and parse */

#if LOG_LEVEL >= 3
        Timer timer(comm);
        timer.start();
#endif

        write_sendbufs(sendbuf_x);

#if LOG_LEVEL >= 3
        timer.stop_and_log("write_sendbufs");
        timer.start();
#endif

        MPI_Wait(&req, MPI_STATUS_IGNORE);

#if LOG_LEVEL >= 3
        timer.stop_and_log("MPI_Wait");
        timer.start();
#endif

        std::swap(sendbuf_x, sendbuf_y);
        std::swap(recvbuf_x, recvbuf_y);

        if (check_complete(recvbuf_x)) {
            status = BATCH_DONE;
            parse_recvbufs(recvbuf_x);
            return;
        }

        MPI_Ialltoall(sendbuf_y, batch_size, MPI_BYTE, recvbuf_y, batch_size, MPI_BYTE, comm, &req);

#if LOG_LEVEL >= 3
        timer.stop_and_log("MPI_Ialltoall");
        timer.start();
#endif

        parse_recvbufs(recvbuf_x);

#if LOG_LEVEL >= 3
        timer.stop_and_log("parse_recvbufs");
#endif

    
    }

    virtual void print_stats(){
        Logger logger(comm);
        logger() << "Round: "<< round << std::endl;
        logger.flush("Print Stats", 0);
    }

    BatchExchanger(MPI_Comm comm, size_t ntasks, size_t batch_size, size_t max_element_size) : 
        comm(comm), batch_size(batch_size), status(BATCH_NOT_INIT), ntasks(ntasks), max_element_size(max_element_size)
    {
        round = 0;
        MPI_Comm_size(comm, &nprocs);
        MPI_Comm_rank(comm, &myrank);
        send_limit = batch_size - sizeof(char) - max_element_size;

    };

    ~BatchExchanger()  
    {
        delete[] sendbuf_x;
        delete[] recvbuf_x;
        delete[] sendbuf_y;
        delete[] recvbuf_y;
    }
};

class LengthExchanger : public BatchExchanger
{
private:
    int nthr_membounded;
    std::vector<std::vector<std::vector<uint32_t>>>& lengths;
    
    std::vector<size_t> current_taskid;
    std::vector<size_t> current_tid;
    std::vector<size_t> current_idx;

    bool write_sendbuf(uint8_t* addr, int procid) override {
        size_t taskid = current_taskid[procid];
        size_t tid = current_tid[procid];
        size_t idx = current_idx[procid];

        size_t cnt = 0;
        while (cnt <= send_limit && taskid < (procid + 1) * ntasks ) {
            size_t n = lengths[tid][taskid].size();
            // need to check if this +1 is good
            size_t n_to_send = std::min(n - idx, (send_limit - cnt) / sizeof(uint32_t) + 1);
            memcpy(addr + cnt, lengths[tid][taskid].data() + idx, n_to_send * sizeof(uint32_t));
            cnt += n_to_send * sizeof(uint32_t);
            idx += n_to_send;

            if (idx == n) {
                tid++;
                idx = 0;
                if(tid == nthr_membounded) {
                    taskid++;
                    tid = 0;
                }
            }
        }

        current_taskid[procid] = taskid;
        current_tid[procid] = tid;
        current_idx[procid] = idx;

        return taskid == (procid + 1) * ntasks;
    }

    std::vector<size_t>& recv_cnt;
    std::vector<std::vector<uint32_t>>& recvbuf;
    std::vector<size_t> recv_idx;
    std::vector<size_t> recv_taskid;
    // Note that here the task id is the "local" task id, which starts from 0 for each process

    void parse_recvbuf(uint8_t* addr, int procid) override {
        size_t taskid = recv_taskid[procid];
        size_t idx = recv_idx[procid];

        if (taskid == ntasks) {
            return;
        }

        size_t cnt = 0;
        while (cnt <= send_limit && taskid < ntasks) {
            size_t n_to_recv = std::min(recv_cnt[taskid + procid * ntasks] - idx, (send_limit - cnt) / sizeof(uint32_t) + 1);
            memcpy(recvbuf[taskid + procid * ntasks].data() + idx, addr + cnt, n_to_recv * sizeof(uint32_t));
            cnt += n_to_recv * sizeof(uint32_t);
            idx += n_to_recv;
            if (idx == recv_cnt[taskid + procid * ntasks]) {
                taskid++;
                idx = 0;
            }
        }

        recv_taskid[procid] = taskid;
        recv_idx[procid] = idx;

    }

public:
    LengthExchanger(MPI_Comm comm, 
                size_t ntasks, 
                size_t batch_size, 
                size_t max_element_size, 
                int nthr_membounded, 
                std::vector<std::vector<std::vector<uint32_t>>>& lengths,
                std::vector<size_t>& recv_cnt,
                std::vector<std::vector<uint32_t>>& recvbuf) : 
        BatchExchanger(comm, ntasks, batch_size, max_element_size), 
        nthr_membounded(nthr_membounded), lengths(lengths), recv_cnt(recv_cnt), recvbuf(recvbuf) {
            current_taskid.resize(nprocs);
            for (int i = 0; i < nprocs; i++) {
                current_taskid[i] = ntasks * i;
            }
            current_tid.resize(nprocs, 0);
            current_idx.resize(nprocs, 0);

            recv_idx.resize(nprocs, 0);
            recv_taskid.resize(nprocs, 0);
            recvbuf.resize(nprocs * ntasks);
            for (int i = 0; i < nprocs * ntasks; i++) {
                recvbuf[i].resize(recv_cnt[i], 0);
            }
        }
};



struct BucketAssistant{
    int ntasks, nprocs;
    KmerSeedBuckets& bucket;
    std::vector<std::vector<size_t>> recv_base;
    std::vector<std::vector<size_t>> current_recv;


    BucketAssistant(int ntasks, int nprocs, KmerSeedBuckets& bucket, std::vector<std::vector<uint32_t>>& lengths) : 
        ntasks(ntasks), nprocs(nprocs), bucket(bucket){
        std::vector<std::vector<size_t>> recv_cnt;

        recv_cnt.resize(nprocs);
        recv_base.resize(nprocs);
        current_recv.resize(nprocs);

        // Turn the lengths into the actual counts
        for(int i = 0; i < nprocs; i++) {
            recv_cnt[i].resize(ntasks, 0);
            for(int j = 0; j < ntasks; j++) {
                recv_cnt[i][j] = std::accumulate(lengths[i * ntasks + j].begin(), lengths[i * ntasks + j].end(), 0) - ( KMER_SIZE - 1 ) * lengths[i * ntasks + j].size();
            }
        }

        recv_base[0].resize(ntasks, 0);
        for(int i = 1; i < nprocs; i++) {
            recv_base[i].resize(ntasks, 0);
            for(int j = 0; j < ntasks; j++) {
                recv_base[i][j] = recv_base[i-1][j] + recv_cnt[i-1][j];
            }
        }

        for(int i = 0; i < nprocs; i++) {
            current_recv[i].resize(ntasks, 0);
        }

        for(int i = 0; i < ntasks; i++) {
            bucket[i].resize(recv_base[nprocs-1][i] + recv_cnt[nprocs-1][i]);
        }
    }    
    
    inline void insert(const int& procid, const int& taskid, const DnaSeq& seq) {
        size_t len = seq.size() - KMER_SIZE + 1;

        auto repmers = TKmer::GetRepKmers(seq);
        size_t base = recv_base[procid][taskid] + current_recv[procid][taskid];

        for (int i = 0; i < len; i++) {
            bucket[taskid][base + i] = KmerSeedStruct(repmers[i]);
        }

        current_recv[procid][taskid] += len;
    }
};

class SupermerExchanger : public BatchExchanger
{
private:
    int nthr_membounded;
    std::vector<std::vector<std::vector<uint32_t>>>& lengths;
    std::vector<std::vector<std::vector<uint8_t>>>& supermers;
    std::vector<size_t> current_taskid;
    std::vector<size_t> current_tid;
    std::vector<size_t> current_idx;
    std::vector<size_t> current_supermer_idx;

    std::vector<size_t> _bytes_sent;

    bool write_sendbuf(uint8_t* addr, int procid) override {
        size_t taskid = current_taskid[procid];
        size_t tid = current_tid[procid];
        size_t idx = current_idx[procid];
        size_t supermer_idx = current_supermer_idx[procid];

        size_t cnt = 0;
        while (cnt <= send_limit && taskid < (procid + 1) * ntasks ) {

            if(idx >= lengths[tid][taskid].size()) {
                tid++;
                idx = 0;
                supermer_idx = 0;
                if(tid == nthr_membounded) {
                    taskid++;
                    tid = 0;
                }
                continue;
            }

            size_t len_bytes = cnt_bytes(lengths[tid][taskid][idx]);
            memcpy(addr + cnt, supermers[tid][taskid].data() + supermer_idx, len_bytes);
            cnt += len_bytes;
            supermer_idx += len_bytes;
            idx++;

        }

        current_taskid[procid] = taskid;
        current_tid[procid] = tid;
        current_idx[procid] = idx;
        current_supermer_idx[procid] = supermer_idx;

        return taskid == (procid + 1) * ntasks;
    }

    std::vector<std::vector<uint32_t>>& recv_length;
    std::vector<size_t> recv_idx;
    std::vector<size_t> recv_taskid;
    std::vector<size_t> recv_cnt;
    BucketAssistant assistant;

    void parse_recvbufs(uint8_t* addr) {
        #pragma omp parallel for num_threads(MAX_THREAD_MEMORY_BOUNDED)
        for (int i = 0; i < nprocs; i++) {
            parse_recvbuf(addr + i * batch_size, i);
        }
    }


    void parse_recvbuf(uint8_t* addr, int procid) override {
        size_t taskid = recv_taskid[procid];
        size_t idx = recv_idx[procid];

        if (taskid == ntasks) {
            return;
        }

        size_t cnt = 0;
        while (cnt <= send_limit && taskid < ntasks) {
            if (idx >= recv_cnt[procid * ntasks + taskid]) {
                taskid++;
                idx = 0;
                continue;
            }
            size_t len = recv_length[procid * ntasks + taskid][idx];
            size_t len_bytes = cnt_bytes(len);

            auto seq = DnaSeq(len, addr+cnt);
            assistant.insert(procid, taskid, seq);

            // memcpy(bucket[procid * ntasks + taskid].data() + idx, addr + cnt, len_bytes);
            cnt += len_bytes;
            idx++;
        }

        recv_taskid[procid] = taskid;
        recv_idx[procid] = idx;

        _bytes_sent[procid] += cnt;

    }


public:
    SupermerExchanger(MPI_Comm comm, 
                size_t ntasks, 
                size_t batch_size, 
                size_t max_element_size, 
                int nthr_membounded, 
                std::vector<std::vector<std::vector<uint32_t>>>& lengths,
                std::vector<std::vector<std::vector<uint8_t>>>& supermers,
                std::vector<size_t>& recv_cnt,
                std::vector<std::vector<uint32_t>>& recv_length,
                KmerSeedBuckets& bucket) : 
        BatchExchanger(comm, ntasks, batch_size, max_element_size), 
        nthr_membounded(nthr_membounded), lengths(lengths), supermers(supermers), recv_cnt(recv_cnt), recv_length(recv_length), 
        assistant(ntasks, nprocs, bucket, recv_length) {
            current_taskid.resize(nprocs);
            for (int i = 0; i < nprocs; i++) {
                current_taskid[i] = ntasks * i;
            }
            current_tid.resize(nprocs, 0);
            current_idx.resize(nprocs, 0);
            current_supermer_idx.resize(nprocs, 0);

            recv_idx.resize(nprocs, 0);
            recv_taskid.resize(nprocs, 0);
            
            _bytes_sent.resize(nprocs, 0);

        }

    void print_stats() override {
        Logger logger(comm);
        
        // Rank 0 Gather the bytes sent by each process
        std::vector<size_t> all_bytes_sent(nprocs * nprocs, 0);
        MPI_Gather(_bytes_sent.data(), nprocs, MPI_UNSIGNED_LONG, all_bytes_sent.data(), nprocs, MPI_UNSIGNED_LONG, 0, comm);

        // Each Rank Calculate its own std deviation
        size_t mean = 0;
        size_t std_dev = 0;
        for(int i = 0; i < nprocs; i++) {
            mean += _bytes_sent[i];
        }
        mean /= nprocs;
        for (int i = 0; i < nprocs; i++) {
            std_dev += (_bytes_sent[i] - mean) * (_bytes_sent[i] - mean);
        }
        std_dev = std::sqrt(std_dev / nprocs);
        logger() << "mean: " << mean << " \t\t std_dev: " << std_dev << " \t\t ratio: " << (double)std_dev / mean  ;
        logger.flush("Sending Stats in Bytes");

        // Calculate the std deviation for all ranks sending
        if(myrank == 0) {
            size_t mean_all = 0;
            size_t std_dev_all = 0;
            for(int i = 0; i < nprocs * nprocs; i++) {
                mean_all += all_bytes_sent[i];
            }
            mean_all /= nprocs * nprocs;
            for (int i = 0; i < nprocs * nprocs; i++) {
                std_dev_all += (all_bytes_sent[i] - mean_all) * (all_bytes_sent[i] - mean_all);
            }
            std_dev_all = std::sqrt(std_dev_all / (nprocs * nprocs));
            logger() << "mean: " << mean_all << " \t\t std_dev: " << std_dev_all << " \t\t ratio: " << (double)std_dev_all / mean_all  ;
        }
        logger.flush("Overall Sending Stats in bytes", 0);

        // Calculate the average length of the supermers
        size_t total_supermer_len = 0;
        size_t total_supermer_cnt = 0;
        for (int i = 0; i < lengths.size(); i++) {
            for (int j = 0; j < lengths[i].size(); j++) {
                total_supermer_len += std::accumulate(lengths[i][j].begin(), lengths[i][j].end(), 0);
                total_supermer_cnt += lengths[i][j].size();
            }
        }

        // Rank 0 gather the total supermer length and count
        size_t all_supermer_len = 0;
        size_t all_supermer_cnt = 0;
        MPI_Reduce(&total_supermer_len, &all_supermer_len, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&total_supermer_cnt, &all_supermer_cnt, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);

        logger() << "Average supermer length: " << (double)all_supermer_len / all_supermer_cnt;
        logger.flush("Supermer Stats", 0);

        logger()<< "Round: "<< round << std::endl;
        logger.flush("Round Stats", 0);
    
    }

};

std::unique_ptr<KmerSeedBuckets> exchange_supermer(ParallelData& data, MPI_Comm comm);


struct KmerParserHandler
{
    int nprocs;
    ReadId readoffset;
    std::vector<std::vector<KmerSeedStruct>>& kmerseeds;

    KmerParserHandler(std::vector<std::vector<KmerSeedStruct>>& kmerseeds, ReadId readoffset) : nprocs(kmerseeds.size()), readoffset(readoffset), kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer, size_t kid, size_t rid)
    {
        kmerseeds[GetKmerOwner(kmer, nprocs)].emplace_back(kmer);
    }
};


template <typename KmerHandler>
void ForeachKmerParallel(const DnaBuffer& myreads, std::vector<KmerHandler>& handlers, int nthreads)
{
    assert(nthreads > 0);
    /*
     * Go through each local read.
     */

    /* cosidering the NUMA effect, we may want the vecs to be NUMA-local*/
    #pragma omp parallel for num_threads(nthreads) 
    for (size_t i = 0; i < myreads.size(); ++i)
    {
        int tid = omp_get_thread_num();
        /*
         * If it is too small then continue to the next one.
         */
        if (myreads[i].size() < KMER_SIZE)
            continue;

        /*
         * Get all the representative k-mer seeds.
         */
        std::vector<TKmer> repmers = TKmer::GetRepKmers(myreads[i]);

        size_t j = 0;

        /*
         * Go through each k-mer seed.
         */
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j)
        {
            handlers[tid](*meritr, j, i);
        }
    }
}

struct Minimizer_Deque {

    /* The first item is the hash value, the second one is the position of the minimizer */
    std::deque<std::pair<uint64_t, int>> deq;

    void remove_minimizer(int pos) {
        while (!deq.empty() && deq.front().second <= pos) {
            deq.pop_front();
        }
    }

    void insert_minimizer(uint64_t hash, int pos) {
        while (!deq.empty() && deq.back().first < hash) {
            deq.pop_back();
        }
        deq.push_back(std::make_pair(hash, pos));
    }

    uint64_t get_current_minimizer() {
        return deq.front().first;
    }
};

int GetMinimizerOwner(const uint64_t& hash, int tot_tasks);

void FindKmerDestinationsParallel(const DnaBuffer& myreads, int nthreads, int tot_tasks, ParallelData& data);

#endif // KMEROPS_H_