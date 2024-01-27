#include<iostream>
#include<mpi.h>
#include<omp.h>

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
        std::cerr << "Compiling Parameters:" << std::endl;
        std::cerr << "      KMER_SIZE: " << KMER_SIZE << std::endl;
        std::cerr << "      THREAD_PER_TASK: " << THREAD_PER_TASK << std::endl;
        std::cerr << "      LOWER_KMER_FREQ: " << LOWER_KMER_FREQ << std::endl;
        std::cerr << "      UPPER_KMER_FREQ: " << UPPER_KMER_FREQ << std::endl;
        std::cerr << "      LOGGING_LEVEL: " << LOG_LEVEL << std::endl;
        std::cerr << "      DEBUG: " << DEBUG << std::endl<< std::endl;
    }

    if (myrank == 0){
        log() << "Runtime Parameters:" << std::endl;
        log() << "      Fasta File: " << std::quoted(fasta_fname)<< std::endl;
        log() << "      Nprocs:" << nprocs << std::endl;
        log() << "      Default Maximum Thread Count Per Process: " << omp_get_max_threads() << std::endl;
    }
    log.flush(log(), 0);

    timer.start();
    FastaIndex index(fasta_fname, MPI_COMM_WORLD);
    ss << "reading " << std::quoted(index.get_faidx_fname()) << " and scattering to all MPI tasks";
    timer.stop_and_log(ss.str().c_str());

    ss.clear(); ss.str("");


    timer.start();
    DnaBuffer mydna = index.getmydna();
    ss << "reading and 2-bit encoding " << std::quoted(index.get_fasta_fname()) << " sequences in parallel";
    timer.stop_and_log(ss.str().c_str());
    ss.clear(); ss.str("");

    /* start kmer counting */
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