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

typedef std::tuple<TKmer, READIDS, POSITIONS, int> KmerListEntry;
typedef std::vector<KmerListEntry> KmerList;

#define DEFAULT_THREAD_PER_TASK 4
#define MAX_THREAD_MEMORY_BOUNDED 8
#define MAX_SEND_BATCH 1000000

struct KmerSeedStruct{
    TKmer kmer;      
    ReadId readid;
    PosInRead posinread;

    KmerSeedStruct(TKmer kmer, ReadId readid, PosInRead posinread) : kmer(kmer), readid(readid), posinread(posinread) {};
    KmerSeedStruct(const KmerSeedStruct& o) : kmer(o.kmer), readid(o.readid), posinread(o.posinread) {};
    KmerSeedStruct(KmerSeedStruct&& o) :  
        kmer(std::move(o.kmer)), readid(std::move(o.readid)), posinread(std::move(o.posinread)) {};
    KmerSeedStruct() {};

    int GetByte(int &i) const { return kmer.getByte(i); }
    bool operator < (const KmerSeedStruct& o) const { return kmer < o.kmer; }
    bool operator == (const KmerSeedStruct& o) const { return kmer == o.kmer; }
    bool operator != (const KmerSeedStruct& o) const { return kmer != o.kmer; }

    KmerSeedStruct& operator=(const KmerSeedStruct& o)
    {
        kmer = o.kmer;
        readid = o.readid;
        posinread = o.posinread;
        return *this;
    }
};

typedef std::vector<std::vector<KmerSeedStruct>> KmerSeedBuckets;
typedef std::vector<std::vector<std::vector<KmerSeedStruct>>> KmerSeedVecs;

std::unique_ptr<KmerSeedBuckets> 
exchange_kmer(const DnaBuffer& myreads,
     MPI_Comm comm,
     int thr_per_task = DEFAULT_THREAD_PER_TASK,
     int max_thr_membounded = MAX_THREAD_MEMORY_BOUNDED);

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = DEFAULT_THREAD_PER_TASK);

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

    size_t batch_sz;        /* batch size in bytes */
    size_t send_limit;      /* the max number of bytes for data */
    size_t ntasks;             /* the number of tasks */
    size_t nthrs;             /* the number of threads */
    size_t* loc_task;     /* the location of the currently sending task cnt for processes */
    size_t* loc_thr;      /* the location of the currently sending thread cnt for processes */
    size_t* loc_idx;      /* the location of the next sending item for processors */

    KmerSeedVecs& kmerseeds;

    /* x is always the one that is currently sending, y is the one that can be used */
    uint8_t* sendtmp_x;
    uint8_t* recvtmp_x;
    uint8_t* sendtmp_y;
    uint8_t* recvtmp_y;

    /* encode and write a kmerseed to the current buffer address */
    void encode_kmerseed_basic(uint8_t* &addr, KmerSeedStruct& kmerseed) {

        TKmer kmer = kmerseed.kmer;
        ReadId readid = kmerseed.readid;
        PosInRead pos = kmerseed.posinread;

        memcpy(addr, kmer.GetBytes(), TKmer::NBYTES);
        memcpy(addr + TKmer::NBYTES, &readid, sizeof(ReadId));
        memcpy(addr + TKmer::NBYTES + sizeof(ReadId), &pos, sizeof(PosInRead));

        addr += (sizeof(TKmer) + sizeof(ReadId) + sizeof(PosInRead));
    }

    /* encoder vars */
    ReadId last_readid;
    PosInRead last_posinread;

    void init_encoder() {
        last_readid = 0;
        last_posinread = 0;
    }

    int ca = 0, cb = 0, cc = 0, cd = 0;
    
    void encode_kmerseed(uint8_t* &addr, KmerSeedStruct& kmerseed) {

        TKmer kmer = kmerseed.kmer;
        ReadId readid = kmerseed.readid;
        PosInRead pos = kmerseed.posinread;

        memcpy(addr, kmer.GetBytes(), TKmer::NBYTES);

        char c = 0;
        if (readid == last_readid) {
            if (pos > last_posinread && pos - last_posinread < (2 << 16))
                c = 3;
            else
                c = 2;
        } else if (readid == last_readid + 1) {
            c = 1;
        };

        switch (c) {
            case 0:
                memcpy(addr + TKmer::NBYTES, &c, sizeof(char));
                memcpy(addr + TKmer::NBYTES + sizeof(char), &pos, sizeof(PosInRead));
                memcpy(addr + TKmer::NBYTES + sizeof(char) + sizeof(PosInRead), &readid, sizeof(ReadId));
                addr += (sizeof(TKmer) + sizeof(char) + sizeof(ReadId) + sizeof(PosInRead));
                ca++;
                break;
            case 1:
                memcpy(addr + TKmer::NBYTES, &c, sizeof(char));
                memcpy(addr + TKmer::NBYTES + sizeof(char), &pos, sizeof(PosInRead));
                addr += (sizeof(TKmer) + sizeof(char) + sizeof(PosInRead));
                cb++;
                break;
            case 2:
                memcpy(addr + TKmer::NBYTES, &c, sizeof(char));
                memcpy(addr + TKmer::NBYTES + sizeof(char), &pos, sizeof(PosInRead));
                addr += (sizeof(TKmer) + sizeof(char) + sizeof(PosInRead));
                cc++;
                break;
            case 3:
                uint16_t pos16 = (uint16_t)(pos - last_posinread);
                memcpy(addr + TKmer::NBYTES, &c, sizeof(char));
                memcpy(addr + TKmer::NBYTES + sizeof(char), &pos16, sizeof(uint16_t));
                addr += (sizeof(TKmer) + sizeof(char) + sizeof(uint16_t));
                cd++;
                break;
        }

        last_readid = readid;
        last_posinread = pos;

    }

    /* 
     * Write sendbuf for a process. returns if it has completed 
     * This function is also responsible for determining the vector to read from
     * Returns if all data has been sent (for this process)
     */
    bool write_sendbuf(uint8_t* addr, int procid) {

        init_encoder();
        size_t task = loc_task[procid];
        size_t thr = loc_thr[procid];
        size_t idx = loc_idx[procid];

        std::vector<KmerSeedStruct> &v = kmerseeds[thr][procid * ntasks + task];
        // std::cout<<kmerseeds[thr].size()<<" "<<procid * ntasks + task<<std::endl;

        if (task >= ntasks) return true;    // all data has been sent. is there a bug in the original version?

        size_t max_idx = v.size();

        uint8_t* addr_limit = addr + send_limit;

        while (addr <= addr_limit && task < ntasks) {
            encode_kmerseed(addr, v[idx]);
            idx++;

            if (idx >= max_idx) {
                idx = 0;
                thr++;
                if (thr >= nthrs) {
                    thr = 0;
                    task++;
                }
                if (task >= ntasks) break;

                v = kmerseeds[thr][procid * ntasks + task];
                max_idx = v.size();
            }
        }

        loc_idx[procid] = idx;
        loc_thr[procid] = thr;
        loc_task[procid] = task;

        return (task >= ntasks);
    }

    /* write the whole sendbuf for one pass. returns if all data has been sent */
    void write_sendbufs(uint8_t* addr) {

        bool finished = true;

        // #pragma omp parallel for num_threads(4)
        for(int i = 0; i < nprocs; i++) {
            bool myfinished = write_sendbuf(addr + batch_sz * i, i);
            
            // #pragma omp critical
            {
                if (!myfinished) finished = false;
            }
        }

        /* write the last byte of each buf which indicates finished or not */
        if (finished) {
            for (int i = 0; i < nprocs; i++)
                addr[batch_sz * i + batch_sz - 1] = 1;
        } else {
            for (int i = 0; i < nprocs; i++)
                addr[batch_sz * i + batch_sz - 1] = 0;
        }

    }


public:

    void print_results(Logger &l) {
        l()<<ca<<" "<<cb<<" "<<cc<<" "<<cd<<std::endl;
        l.flush("compress results");
    }
    BatchSender(MPI_Comm comm, size_t ntasks, size_t nthrs, size_t batch_size, KmerSeedVecs& kmerseeds) : 
        comm(comm), batch_sz(batch_size), kmerseeds(kmerseeds), status(BATCH_NOT_INIT), ntasks(ntasks), nthrs(nthrs)
    {
        /*
         * Suppose buffer has a max size of n
         * n - size(kmer) - size(readid) - size(posinread) - sizeof(char) = max start location
         * which means that when the pointer is leq this value, we must continue to store items
         * if the pointer is greater than the value, we cannot store items anymore
         * the last byte is reserved for the flag, indicating if a process has sent all its data .
         * 
         * We also need to take care: what if I have sent all my data but some other processes haven't?
         * Well, we just sent random bits to the process. The Important thing is that the last flag is always set.
         * The dst process know how many kmers we're sending, so it will limit the number of kmers it reads.
         */

        send_limit = batch_sz - sizeof(char) - TKmer::NBYTES - sizeof(ReadId) - sizeof(PosInRead) - sizeof(char);

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

    void init()
    {

        status = BATCH_SENDING;
        write_sendbufs(sendtmp_x);
        MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);

    }

    uint8_t* progress() {
        if (status != BATCH_SENDING) return 0;

        #if LOG_LEVEL >= 3
        MPITimer timer(comm);
        timer.start();
        #endif

        write_sendbufs(sendtmp_y);

        #if LOG_LEVEL >= 3
        timer.stop_and_log("write_sendbufs");
        timer.start();
        #endif

        MPI_Wait(&req, MPI_STATUS_IGNORE);

        #if LOG_LEVEL >= 3
        timer.stop_and_log("MPI_Wait");
        #endif

        std::swap(sendtmp_x, sendtmp_y);
        std::swap(recvtmp_x, recvtmp_y); 

        // see if completed
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


        MPI_Ialltoall(sendtmp_x, batch_sz, MPI_BYTE, recvtmp_x, batch_sz, MPI_BYTE, comm, &req);
        return recvtmp_y;
    }
};


class BatchStorer
{
    KmerSeedBuckets& kmerseeds;
    int taskcnt, nprocs;
    size_t batch_size;
    size_t send_limit;
    size_t* expected_recvcnt;
    size_t* recvcnt;
    int* recv_task_id;

    /* it should be inlined by default */
    inline void decode_basic(uint8_t* &addr, std::vector<KmerSeedStruct> &kmerseeds) {
        TKmer kmer((void*)addr);                
        ReadId readid = *((ReadId*)(addr + TKmer::NBYTES));
        PosInRead pos = *((PosInRead*)(addr + TKmer::NBYTES + sizeof(ReadId)));
        kmerseeds.emplace_back(kmer, readid, pos);
        addr += TKmer::NBYTES + sizeof(ReadId) + sizeof(PosInRead);
    }

    /* decoder vars */
    ReadId last_readid;
    PosInRead last_posinread;

    void init_decoder() {
        last_readid = 0;
        last_posinread = 0;
    }

    inline void decode(uint8_t* &addr, std::vector<KmerSeedStruct> &kmerseeds) {
        TKmer kmer((void*)addr);
        
        char c = *((char*)(addr + TKmer::NBYTES));
        addr = addr + TKmer::NBYTES + sizeof(char);

        ReadId readid = 0;
        PosInRead pos = 0;

        switch (c) {
            case 0:
                pos = (*((PosInRead*)(addr)));
                readid = *((ReadId*)(addr + sizeof(PosInRead)));
                addr += (sizeof(ReadId) + sizeof(PosInRead));
                break;
            case 1:
                pos = (*((PosInRead*)(addr)));
                readid = last_readid + 1;
                addr += (sizeof(PosInRead));
                break;
            case 2:
                pos = (*((PosInRead*)(addr)));
                readid = last_readid;
                addr += (sizeof(PosInRead));
                break;
            case 3:
                uint16_t pos16 = *((uint16_t*)(addr));
                pos = last_posinread + pos16;
                readid = last_readid;
                addr += (sizeof(uint16_t));
                break;
        }

        kmerseeds.emplace_back(kmer, readid, pos);

        last_readid = readid;
        last_posinread = pos;
    }

public:
    BatchStorer(KmerSeedBuckets& kmerseeds,  int taskcnt, int nprocs, size_t* expected_recvcnt, size_t batch_size) : 
        kmerseeds(kmerseeds), taskcnt(taskcnt), nprocs(nprocs), expected_recvcnt(expected_recvcnt), batch_size(batch_size)
        {
            send_limit = batch_size - sizeof(char) - TKmer::NBYTES - sizeof(ReadId) - sizeof(PosInRead) - sizeof(char);
            recvcnt = new size_t[nprocs];
            recv_task_id = new int[nprocs];
            memset(recvcnt, 0, sizeof(size_t) * nprocs);
            memset(recv_task_id, 0, sizeof(int) * nprocs);
        }

    ~BatchStorer()
    {
        
        delete[] recvcnt;
        delete[] recv_task_id;
        
    }

    void store(uint8_t* recvbuf)
    {
        #if LOG_LEVEL >= 3
        MPITimer timer(MPI_COMM_WORLD);
        timer.start();
        #endif

        for(int i = 0; i < nprocs; i++) 
        {
            init_decoder();
            if (recv_task_id[i] >= taskcnt) continue;      // all data for this task has been received. all the rest is padding
            int working_task = recv_task_id[i];
            size_t cnt = recvcnt[i];
            size_t max_cnt = expected_recvcnt[ taskcnt * i + working_task ];
            uint8_t* addrs2read = recvbuf + batch_size * i;
            uint8_t* addr_limit = addrs2read + send_limit;

            while(addrs2read <= addr_limit && working_task < taskcnt) {
                decode(addrs2read, kmerseeds[working_task]);
                cnt++;

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
        timer.stop_and_log("store");
        #endif
    }
};

struct KmerParserHandler
{
    int nprocs;
    ReadId readoffset;
    std::vector<std::vector<KmerSeedStruct>>& kmerseeds;

    KmerParserHandler(std::vector<std::vector<KmerSeedStruct>>& kmerseeds, ReadId readoffset) : nprocs(kmerseeds.size()), readoffset(readoffset), kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer, size_t kid, size_t rid)
    {
        kmerseeds[GetKmerOwner(kmer, nprocs)].emplace_back(kmer, static_cast<ReadId>(rid) + readoffset, static_cast<PosInRead>(kid));
    }
};


template <typename KmerHandler>
void ForeachKmerParallel(const DnaBuffer& myreads, std::vector<KmerHandler>& handlers, int nthreads)
{
    assert(nthreads > 0);
    /*
     * Go through each local read.
     */

    // yfli: TODO: test the performance of different number of threads. This might be a memory bandwidth issue.

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