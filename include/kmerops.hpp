#ifndef KMEROPS_H_
#define KMEROPS_H_

#include <mpi.h>
#include <omp.h>
#include <deque>
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

#define MAX_THREAD_MEMORY_BOUNDED 8 

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

void
prepare_supermer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int thr_per_task = THREAD_PER_TASK,
     int max_thr_membounded = MAX_THREAD_MEMORY_BOUNDED);

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = THREAD_PER_TASK);

int GetKmerOwner(const TKmer& kmer, int nprocs);

void print_kmer_histogram(const KmerList& kmerlist, MPI_Comm comm);

class ParallelData{
    private:
        int nprocs;
        int ntasks;
        int nthr_membounded;

        std::vector<std::vector<std::vector<int>>> destinations;
        std::vector<std::vector<std::vector<uint32_t>>> lengths;
        std::vector<std::vector<std::vector<bool>>> supermers;
        std::vector<std::vector<uint64_t>> readids;

    public:
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

        std::vector<std::vector<bool>>& get_my_supermers(int tid) {
            return supermers[tid];
        }

        std::vector<uint64_t>& get_my_readids(int tid) {
            return readids[tid];
        }
};

inline int pad_bytes(const int& len) {
    return (4 - len % 4) * 2;
}

struct SupermerEncoder{
    std::vector<std::vector<uint32_t>>& lengths;
    std::vector<std::vector<bool>>& supermers;
    int max_supermer_len;

    SupermerEncoder(std::vector<std::vector<uint32_t>>& lengths, 
            std::vector<std::vector<bool>>& supermers, 
            int max_supermer_len) : 
            lengths(lengths), supermers(supermers), max_supermer_len(max_supermer_len) {};


    /* std::vector<bool> is a special container that takes 1 bit per element instead of 1 byte */
    void copy_bits(std::vector<bool>& dst, const uint8_t* src, uint64_t start_pos, int len){

        dst.reserve(len*2);

        for (int i = start_pos; i < start_pos + len; i++) {
            /* the order is confusing, just make sure we keep it consistent*/
            bool first_bit = (src[i/4] >> (7 - 2*(i%4))) & 1;
            bool second_bit = (src[i/4] >> (6 - 2*(i%4))) & 1;
            dst.push_back(first_bit);
            dst.push_back(second_bit);
        }

        /* pad the last byte */
        int pad = pad_bytes(len);
        for (int i = 0; i < pad; i++) {
            dst.push_back(0);
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
                last_dst = dest[i];
                cnt = 0;
                start_pos = i;
            }

            /* increment the counter */
            cnt++;
        }
    }
};

class BatchSender
{
public:
    enum Status
    {
        BATCH_NOT_INIT = 0,
        BATCH_SENDING = 1,
        BATCH_DONE = 2
    };

private:
    MPI_Comm comm;
    int nprocs;
    int myrank;

    Status status;
    MPI_Request req;

    size_t batch_sz;            /* maxium batch size in bytes */
    size_t send_limit;          /* maxium size of meaningful data in bytes */
    size_t ntasks;              /* number of tasks */
    size_t nthrs;               /* number of threads */
    size_t* loc_task;           /* the location of the currently sending (destination) task id for processes */
    size_t* loc_thr;            /* the location of the currently sending (local) thread id for processes */
    size_t* loc_idx;            /* the location of the next sending item for processors */

    KmerSeedVecs& kmerseeds;    /* the kmers to be sent */

    /* 
     * two sending / receiving buffers are used to overlap computation and communication
     * x is always the buffer that is currently sending, y is the one that can be overwritten 
     */
    uint8_t* sendtmp_x;
    uint8_t* recvtmp_x;
    uint8_t* sendtmp_y;
    uint8_t* recvtmp_y;

    /* encode and write a kmerseed to the current buffer address */
    void encode_kmerseed_basic(uint8_t* &addr, KmerSeedStruct& kmerseed) {
        TKmer kmer = kmerseed.kmer;
        memcpy(addr, kmer.GetBytes(), TKmer::NBYTES);
        addr += (sizeof(TKmer));
    }

    /* 
     * Write sendbuf for one destination process. 
     * Returns true if all kmers to this process have been sent.
     */
    bool write_sendbuf(uint8_t* addr, int procid);

    /* 
     * Write sendbuf for all destination processes
     */
    void write_sendbufs(uint8_t* addr);


public:

    /* we may want to log something */
    void print_results(Logger &l) {}

    BatchSender(MPI_Comm comm, size_t ntasks, size_t nthrs, size_t batch_size, KmerSeedVecs& kmerseeds) : 
        comm(comm), batch_sz(batch_size), kmerseeds(kmerseeds), status(BATCH_NOT_INIT), ntasks(ntasks), nthrs(nthrs)
    {
        /*
         * Suppose buffer has a max size of n
         * n - size(kmer) - size(readid) - size(posinread) - sizeof(char) = max start location
         * which means that when the pointer is leq this value (and we still have kmers to send), we must continue to encode and write
         * Or, if the pointer is greater than the value, we cannot store items anymore
         * the last byte is reserved for the flag, indicating if the current process has sent all its data .
         * 
         * What if I have sent all my data but some other processes haven't?
         * Well, just send random bits. The important thing is that the flag is always set.
         * The destination process knows how many kmers we're sending, so it won't decode meaningless data.
         */

        send_limit = batch_sz - sizeof(char) - TKmer::NBYTES - sizeof(char);

        MPI_Comm_size(comm, &nprocs);
        MPI_Comm_rank(comm, &myrank);

        sendtmp_x = new uint8_t[batch_sz * nprocs];
        recvtmp_x = new uint8_t[batch_sz * nprocs];

        sendtmp_y = new uint8_t[batch_sz * nprocs];
        recvtmp_y = new uint8_t[batch_sz * nprocs];

        loc_task = new size_t[nprocs];
        loc_thr = new size_t[nprocs];
        loc_idx = new size_t[nprocs];
        memset(loc_task, 0, sizeof(size_t) * nprocs);
        memset(loc_thr, 0, sizeof(size_t) * nprocs);
        memset(loc_idx, 0, sizeof(size_t) * nprocs);

    };

    ~BatchSender()  
    {
        if (sendtmp_x != nullptr) delete[] sendtmp_x;
        if (recvtmp_y != nullptr) delete[] recvtmp_x;
        if (sendtmp_y != nullptr) delete[] sendtmp_y;
        if (recvtmp_y != nullptr) delete[] recvtmp_y;
        if (loc_task != nullptr) delete[] loc_task;
        if (loc_thr != nullptr) delete[] loc_thr;
        if (loc_idx != nullptr) delete[] loc_idx;

    }

    Status get_status() { return status; }

    /* initialize the sender, send the first batch of data */
    void init() {
        status = BATCH_SENDING;
        write_sendbufs(sendtmp_x);
        MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);

    }

    /* 
     * progress to another round 
     * returns the receiving buffer from the last round
     */
    uint8_t* progress();
    uint8_t* progress_basic_overlap();
    uint8_t* progress_basic_a();
    void progress_basic_b();
};


class BatchReceiver
{
    KmerSeedBuckets& kmerseeds;
    int taskcnt, nprocs;
    size_t batch_size;
    size_t send_limit;
    size_t* expected_recvcnt;           /* the expected total recv counts */
    size_t* recvcnt;                    /* current recv counts */
    int* recv_task_id;                  /* current receiving task ids*/

    /* decode the raw bytes into kmer structs */
    inline void decode_basic(uint8_t* &addr, std::vector<KmerSeedStruct> &kmerseeds) {
        TKmer kmer((void*)addr);                
        kmerseeds.emplace_back(kmer);
        addr += TKmer::NBYTES;
    }


public:
    BatchReceiver(KmerSeedBuckets& kmerseeds,  int taskcnt, int nprocs, size_t* expected_recvcnt, size_t batch_size) : 
        kmerseeds(kmerseeds), taskcnt(taskcnt), nprocs(nprocs), expected_recvcnt(expected_recvcnt), batch_size(batch_size)
        {
            send_limit = batch_size - sizeof(char) - TKmer::NBYTES - sizeof(char);
            recvcnt = new size_t[nprocs];
            recv_task_id = new int[nprocs];
            memset(recvcnt, 0, sizeof(size_t) * nprocs);
            memset(recv_task_id, 0, sizeof(int) * nprocs);
        }

    ~BatchReceiver()
    {
        
        delete[] recvcnt;
        delete[] recv_task_id;
        
    }

    /*
     * Receive data from the sender.
     * Partition them into buckets
     */
    void receive(uint8_t* recvbuf);
};

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