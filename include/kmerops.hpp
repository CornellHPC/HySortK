#ifndef KMEROPS_H_
#define KMEROPS_H_

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

/* thread count used for sending and receiving for each upcxx process */
#define THREADCNT_SEND 4
#define THREADCNT_RECV 1
/* send buffer in bytes */
#define BATCH_SIZE 100000
/* thread count used for (openmp) radix sort. this will determine the task count */
#define DEFAULT_THREAD_PER_TASK 4

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

void store_into_local(const int& tid, const size_t& buf_sz, 
    const upcxx::view<uint8_t>& buf_in_rpc, const bool last_send);

std::unique_ptr<KmerSeedBuckets> exchange_kmer(const DnaBuffer& myreads,
     int thr_per_task = DEFAULT_THREAD_PER_TASK);

std::unique_ptr<KmerList> filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = DEFAULT_THREAD_PER_TASK);

int GetKmerOwner(const TKmer& kmer, int nprocs);

void print_kmer_histogram(const KmerList& kmerlist);


class KmerParserHandler
{
    int nprocs;
    int ntasks;

    size_t batch_sz;
    ReadId readoffset;

    std::vector<std::vector<uint8_t>> sendbufs;
    std::vector<size_t> buf_sz;
    upcxx::future<> futures;

    int cnt;
    double time_tot;
    
    void encode(std::vector<uint8_t>& t, const TKmer& kmer, PosInRead posinread, ReadId readid, size_t& buf_sz);
    upcxx::future<> send(int dst, bool last_send = false);

public:
    KmerParserHandler(int ntasks, int nprocs, size_t batch_sz, ReadId readoffset) : 
        ntasks(ntasks), nprocs(nprocs), readoffset(readoffset), batch_sz(batch_sz) {
            buf_sz.resize(nprocs * ntasks);
            sendbufs.resize(nprocs * ntasks);
            for (std::vector<uint8_t>& buf : sendbufs) {
                buf.resize(batch_sz + sizeof(TKmer) + sizeof(ReadId) + sizeof(PosInRead));
            }
            futures = upcxx::make_future();
        }


    /* send all the remaining in the buffer and wait for all sends to complete */
    void wait();
    void operator()(const TKmer& kmer, size_t kid, size_t rid);
    void data() {
        std::cout << "total send count: " << cnt << std::endl;
        std::cout << "avg time for calling rpc" << time_tot / cnt<< std::endl;
    }
};

template <typename KmerHandler>
void ForeachKmer(const DnaBuffer& myreads, KmerHandler& handler, size_t st, size_t ed) {
    size_t i;

    for (i = st; i < ed; ++i) {
        if (myreads[i].size() < KMER_SIZE)
            continue;

        std::vector<TKmer> repmers = TKmer::GetRepKmers(myreads[i]);
        size_t j = 0;
        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j) {
            handler(*meritr, j, i);
        }
    }
}

#endif // KMEROPS_H_