#include<iostream>
#include<mpi.h>

#include "logger.hpp"
#include "timer.hpp"
#include "fastaindex.hpp"
#include "dnabuffer.hpp"
#include "dnaseq.hpp"
#include "kmerops.hpp"
#include "compiletime.h"

std::string fasta_fname;

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    Logger log(MPI_COMM_WORLD);
    Timer timer(MPI_COMM_WORLD);
    std::ostringstream ss;

    if (argc < 2){
        std::cerr << "Usage: " << argv[0] << " <fasta file>" << std::endl;
        exit(1);
    }

    fasta_fname = argv[1];

    int myrank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myrank == 0){
        log() << "Fasta File: " << std::quoted(fasta_fname)<< std::endl;
        log() << "Kmer Size: "<< KMER_SIZE << std::endl;
        log() << "nprocs:" << nprocs << std::endl;

    }
    log.flush(log(), 0);

    timer.start();
    FastaIndex index(fasta_fname, MPI_COMM_WORLD);
    ss << "reading " << std::quoted(index.get_faidx_fname()) << " and scattering to all MPI tasks";
    timer.stop_and_log(ss.str().c_str());

    ss.clear(); ss.str("");

    /*
        * DnaBuffer @mydna stores the compressed read sequences assigned
        * to it by the .fai index file, as determined by @index.
        */
    timer.start();
    DnaBuffer mydna = index.getmydna();
    ss << "reading and 2-bit encoding " << std::quoted(index.get_fasta_fname()) << " sequences in parallel";
    timer.stop_and_log(ss.str().c_str());
    ss.clear(); ss.str("");

    timer.start();
    auto bucket = exchange_kmer(mydna, MPI_COMM_WORLD);
    timer.stop_and_log("exchange_kmer");

    timer.start();
    auto kmerlist = filter_kmer(bucket);
    timer.stop_and_log("filter_kmer");

    print_kmer_histogram(*kmerlist, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}