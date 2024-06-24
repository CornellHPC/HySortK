#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include <mpi.h>
#include "dnabuffer.hpp"
#include "compiletime.h"

namespace hysortk {

class FastaIndex
{
public:
    typedef struct { size_t len, pos, bases; } Record;

    FastaIndex(const std::string& fasta_fname, MPI_Comm comm);

    MPI_Comm getcomm() const { return comm; }
    std::string get_fasta_fname() const { return fasta_fname; }
    std::string get_faidx_fname() const { return fasta_fname + ".fai"; }

    size_t gettotrecords() const { return readdispls.back(); }
    size_t getreadcount(size_t i) const { return static_cast<size_t>(readcounts[i]); }
    size_t getreaddispl(size_t i) const { return static_cast<size_t>(readdispls[i]); }
    size_t getmyreadcount() const { int t; MPI_Comm_rank(comm, &t); return getreadcount(t); }
    size_t getmyreaddispl() const { int t; MPI_Comm_size(comm, &t); return getreaddispl(t); }
    int getreadowner(size_t i) const;

    std::vector<size_t> getmyreadlens() const;

    const std::vector<Record>& getmyrecords() const { return myrecords; }
    const std::vector<MPI_Count_t> getreadcounts() const { return readcounts; }
    const std::vector<MPI_Offset_t> getreaddispls() const { return readdispls; }

    DnaBuffer getmydna() const;
    void log(const DnaBuffer& buffer) const;

    static Record get_faidx_record(const std::string& line, std::string& name);

    std::vector<std::string> bcastnames();

private:
    MPI_Comm comm;
    std::vector<Record> myrecords; /* records for the reads local processor is responsible for */
    std::vector<Record> rootrecords; /* records for all the reads across all processors (parsed by root rank in constructor) */
    std::vector<MPI_Count_t> readcounts; /* number of reads assigned to each processor. Each processor gets a copy. |readcounts| == nprocs */
    std::vector<MPI_Offset_t> readdispls; /* displacement counts for reads across all processors. Each processor gets a copy. |readdispls| == nprocs+1 */
    std::string fasta_fname; /* FASTA file name */
    std::vector<std::string> rootnames;
    int nprocs, myrank;

    void getpartition(std::vector<MPI_Count_t>& sendcounts);
};

} // namespace hysortk

#endif