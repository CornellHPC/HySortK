#ifndef KMEROPS_H_
#define KMEROPS_H_

#include <mpi.h>
#include <omp.h>
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

std::unique_ptr<KmerSeedBuckets> 
exchange_kmer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int thr_per_task = THREAD_PER_TASK,
     int max_thr_membounded = MAX_THREAD_MEMORY_BOUNDED);

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = THREAD_PER_TASK);

int GetKmerOwner(const TKmer& kmer, int nprocs);

void print_kmer_histogram(const KmerList& kmerlist, MPI_Comm comm);


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

#endif // KMEROPS_H_