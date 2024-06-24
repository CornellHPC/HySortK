#ifndef HYSORTK_H_
#define HYSORTK_H_

#include "dnabuffer.hpp"
#include "kmer.hpp"
#include <mpi.h>

namespace hysortk {

std::shared_ptr<DnaBuffer> read_dna_buffer(const std::string& fasta_fname, MPI_Comm comm);

std::unique_ptr<KmerListS> kmer_count(const DnaBuffer& mydna, MPI_Comm comm);

void print_kmer_histogram(const KmerListS& kmerlist, MPI_Comm comm);

void write_output_file(const KmerListS& kmerlist, const std::string& output_dir, MPI_Comm comm);

} // namespace hysortk


#endif
