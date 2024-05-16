#include "fastaindex.hpp"
#include "logger.hpp"
#include "timer.hpp"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <cassert>

using Record = typename FastaIndex::Record;

Record FastaIndex::get_faidx_record(const std::string& line, std::string& name)
{
    /*
     * Read a line from a FASTA index file into a record object.
     */
    Record record;
    std::istringstream(line) >> name >> record.len >> record.pos >> record.bases;
    return record;
}

int FastaIndex::getreadowner(size_t i) const
{
    /*
     * We want to find the first processor rank that owns the read with id @i.
     * For each rank i in [0..nprocs-1], readdispls[i] is the number of reads owned across
     * processors [0..i), or the same: the global id of the first read owned by processor i.
     * readdispls[nprocs] is defined as the total number of reads in the FASTA.
     *
     * std:upper_bound(readdispls.cbegin(), readdispls.cend(), i) returns an
     * iterator pointing to the first element in the range [0..nprocs] such that
     * @i < element. This iterator is therefore pointing to the entry of the first
     * processor rank holding a read id greater than @i. It follows that the
     * processor rank just before it owns the read @i.
     */

    auto iditr = std::upper_bound(readdispls.cbegin(), readdispls.cend(), static_cast<uint64_t>(i));

    iditr--;

    return static_cast<int>(iditr - readdispls.cbegin());
}

void FastaIndex::getpartition(std::vector<MPI_Count_t>& sendcounts)
{
    assert(sendcounts.size() == (uint64_t)nprocs);

    size_t totbases = std::accumulate(rootrecords.begin(), rootrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    size_t numreads = rootrecords.size();
    double avgbasesperproc = static_cast<double>(totbases) / nprocs;

    /*
     * Coming up with the optimal partitioning of sequences weighted by their length
     * is NP-hard, and I can't imagine that it is NP-complete also (how would you verify
     * in polynomial time?). Tried looking for the name of this optimization problem
     * but couldn't find (didn't look too hard because its not too important for this).
     * It's basically a variation on multiway number partition, except where the divisions
     * must be ordered. The following seems like a fine approxmimation. Note that the
     * the last processor tends to get more than the average amount of data.
     */

    size_t readid = 0;

    for (int i = 0; i < nprocs-1; ++i)
    {
        size_t basessofar = 0;
        size_t startid = readid;

        /*
         * Keep going through reads until the next one puts us over
         * the average bases per processor.
         */
        do
        {
            basessofar += rootrecords[readid].len;
            readid++;
        } while (readid < numreads && basessofar + rootrecords[readid].len < avgbasesperproc);



        size_t readssofar = readid - startid;
        assert(readssofar >= 1); /* TODO: come up with a recovery strategy here */
        assert(readid <= numreads); /* TODO: come up with a recovery strategy here */

        sendcounts[i] = readssofar;
    }

    /*
     * Last processor gets the remaining reads.
     */
    sendcounts.back() = numreads - readid;
}

FastaIndex::FastaIndex(const std::string& fasta_fname, MPI_Comm comm) :  fasta_fname(fasta_fname), comm(comm)
{

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);
    readcounts.resize(nprocs);



    /*
     * Root processor responsible for reading and parsing FASTA
     * index file "{fasta_fname}.fai" into one record per sequence.
     */
    if (myrank == 0)
    {
        std::string line, name;
        std::ifstream filestream(get_faidx_fname());

        while (std::getline(filestream, line))
        {
            rootrecords.push_back(get_faidx_record(line, name));
            rootnames.push_back(name);
        }

        filestream.close();

        /*
         * Compute load-balanced read partitioning on the root processor.
         */
        getpartition(readcounts);
    }

    /*
     * All processors get a copy of the read counts.
     */
    MPI_Bcast(readcounts.data(), nprocs, MPI_COUNT_TYPE, 0, comm);

    /*
     * And the read displacements.
     */
    readdispls.resize(nprocs);
    std::exclusive_scan(readcounts.begin(), readcounts.end(), readdispls.begin(), static_cast<MPI_Count_t>(0));

    /*
     * It is useful for the displacements to store the total number
     * of reads in the last position.
     */
    // std::cout<<"Counts:"<<readcounts.back()<<std::endl;
    // std::cout<<"Displacements:"<<readdispls.back()<<std::endl;

    readdispls.push_back(readdispls.back() + readcounts.back());

    /*
     * To prevent confusion for the reader: the broadcasting of readcounts
     * to every processor and the parallel "re"-computation of readdispls
     * on each processor is not necessary for performing the scatter operation
     * below, however we still do it because every processor will need to know
     * those things later.
     */
    myrecords.resize(readcounts[myrank]);

    /*
     * Each record is represented with three numbers (read length,
     * FASTA position, FASTA line width), so we create an MPI datatype
     * to communicate each record as a single unit.
     */
    MPI_Datatype faidx_dtype_t;
    MPI_Type_contiguous(3, MPI_UNSIGNED_LONG_LONG, &faidx_dtype_t); // Is this correct?
    MPI_Type_commit(&faidx_dtype_t);

    /*
     * Scatter the records according to the load-balanced read partitioning.
     */

    MPI_Scatterv(rootrecords.data(), readcounts.data(), readdispls.data(), faidx_dtype_t, myrecords.data(), readcounts[myrank], faidx_dtype_t, 0, comm);

    MPI_Type_free(&faidx_dtype_t);


    #if LOG_LEVEL >= 2
    Logger logger(comm);
    size_t mytotbases = std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    size_t totbases;
    MPI_Allreduce(&mytotbases, &totbases, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    double percent_proportion = (static_cast<double>(mytotbases) / totbases) * 100.0;
    logger() << " is responsible for sequences " << Logger::readrangestr(readdispls[myrank], readcounts[myrank]) << " (" << mytotbases << " nucleotides, " << std::fixed << std::setprecision(3) << percent_proportion << "%)";
    logger.flush("Fasta index construction:");
    #endif
}

std::vector<size_t> FastaIndex::getmyreadlens() const
{
    /*
     *  Because we store sequence information using "array of structs"
     *  instead of "structs of arrays", if we want a linear array
     *  of read lengths we have to unpack them.
     */
    std::vector<size_t> readlens(getmyreadcount());
    std::transform(myrecords.cbegin(), myrecords.cend(), readlens.begin(), [](const auto& record) { return record.len; });
    return readlens;
}

DnaBuffer FastaIndex::getmydna() const
{
    /*
     * Allocate local sequence buffer.
     */
    auto readlens = getmyreadlens(); /* vector of local read lengths */
    size_t bufsize = DnaBuffer::computebufsize(readlens); /* minimum number of bytes needed to 2-bit encode all the local reads */
    DnaBuffer dnabuf(bufsize); /* initialize dnabuf by allocating @bufsize bytes */
    size_t numreads = readlens.size(); /* number of local reads */

    MPI_Offset startpos; /* the FASTA position that starts my local chunk of reads */
    MPI_Offset endpos; /* the FASTA position that ends my local chunk of reads (exclusive) */
    MPI_Offset filesize; /* the total size of the FASTA */
    MPI_Offset readbufsize; /* endpos - startpos */
    MPI_File fh;

    /*
     * FASTA file will be read using MPI collective I/O.
     */
    MPI_File_open(comm, get_fasta_fname().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &filesize);

    /*
     * Get start and end coordinates within FASTA of the sequences
     * this processor requires.
     */
    startpos = myrecords.front().pos;
    endpos = myrecords.back().pos + myrecords.back().len + (myrecords.back().len / myrecords.back().bases);
    if (endpos > filesize) endpos = filesize;

    /*
     * Allocate a char buffer to read my FASTA chunk into to,
     * and then do the reading.
     */
    readbufsize = endpos - startpos;
    std::unique_ptr<char[]> readbuf(new char[readbufsize]);

    /*
     * Every processor that calls @getmydna reads in its assigned chunk of data
     * into its temporary local storage buffer (meant for raw contents of file).
     */
    // MPI_File_read_at_all(fh, startpos, &readbuf[0], readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    // MPI_File_close(&fh);

    std::ifstream file(get_fasta_fname().c_str(), std::ios::binary);
    file.seekg(startpos, std::ios::beg);
    file.read(&readbuf[0], readbufsize);
    // MPI_FILE_READ_AT_ALL(fh, startpos, &readbuf[0], readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    file.close();

    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    size_t maxlen = *std::max_element(readlens.begin(), readlens.end());

    /*
     * ASCII sequences are first read into this temporary char buffer
     * before they are compressed into the sequence buffer.
     */
    std::unique_ptr<char[]> tmpbuf(new char[maxlen]);

    MPI_Barrier(comm);
    double elapsed = -MPI_Wtime();

    /*
     * Go through each local FASTA record.
     */
    for (auto itr = myrecords.cbegin(); itr != myrecords.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = &tmpbuf[0];

        /*
         * Read ASCII FASTA sequence into the temoprary buffer.
         */
        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        /*
         * DnaBuffer automatically 2-bit encodes the ASCII sequence
         * and pushes it onto its local stack of sequences.
         */
        dnabuf.push_back(&tmpbuf[0], itr->len);
    }

    elapsed += MPI_Wtime();

    #if LOG_LEVEL >= 2
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(comm);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.flush("FASTA parsing rates (DnaBuffer):");
    #endif

    return dnabuf;
}
