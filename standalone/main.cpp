#include <iostream>
#include<fstream>
#include <mpi.h>
#include <omp.h>
#include "hysortk.hpp"
#include "logger.hpp"
#include "timer.hpp"

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    hysortk::Logger log(MPI_COMM_WORLD);
    hysortk::Timer timer(MPI_COMM_WORLD);
    std::ostringstream ss;

    if (argc < 2){
        std::cerr << "Usage: " << argv[0] << " <fasta file> <output dir>(Optional)" << std::endl;
        exit(1);
    }

    std::string fasta_fname;
    std::string output_dir = "";

    fasta_fname = argv[1];

    if (argc == 3){
        output_dir = argv[2];
    }

    int myrank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myrank == 0){
        log() << "Compiling Parameters:" << std::endl;
        log() << "      KMER_SIZE: " << KMER_SIZE << std::endl;
        log() << "      EXTENSION: " << EXTENSION << std::endl;
        log() << "      MINIMIZER_SIZE: " << MINIMIZER_SIZE << std::endl;
        log() << "      THREAD_PER_WORKER: " << THREAD_PER_WORKER << std::endl;
        log() << "      MAX_THREAD_MEMORY_BOUNDED: " << MAX_THREAD_MEMORY_BOUNDED << std::endl;
        log() << "      LOWER_KMER_FREQ: " << LOWER_KMER_FREQ << std::endl;
        log() << "      UPPER_KMER_FREQ: " << UPPER_KMER_FREQ << std::endl;
        log() << "      LOGGING_LEVEL: " << LOG_LEVEL << std::endl;
        log() << "      DEBUG: " << DEBUG << std::endl;
        log() << "      MAX_SEND_BATCH: " << MAX_SEND_BATCH << std::endl;
        log() << "      AVG_TASK_PER_WORKER: " << AVG_TASK_PER_WORKER << std::endl;
        log() << "      SORT (0: runtime decision, 1: PARADIS, 2: RADULS): " << SORT << std::endl;
        log() << "      DISPATCH_UPPER_COE: " << DISPATCH_UPPER_COE << std::endl;
        log() << "      DISPATCH_STEP: " << DISPATCH_STEP << std::endl;
        log() << "      UNBALANCED_RATIO: " << UNBALANCED_RATIO << std::endl << std::endl;

        log() << "Runtime Parameters:" << std::endl;
        log() << "      Fasta File: " << std::quoted(fasta_fname)<< std::endl;
        log() << "      Output Directory: " << std::quoted(output_dir) << std::endl;
        log() << "      Nprocs:" << nprocs << std::endl;
        log() << "      Default Maximum Thread Count Per Process: " << omp_get_max_threads() << std::endl;
    }
    log.flush(log(), 0);

    auto dna = hysortk::read_dna_buffer(fasta_fname, MPI_COMM_WORLD);

    auto kmer_list = hysortk::kmer_count(dna, MPI_COMM_WORLD);

    hysortk::print_kmer_histogram(*kmer_list, MPI_COMM_WORLD);

    if (argc == 3){
        hysortk::write_output_file(*kmer_list, output_dir, MPI_COMM_WORLD);
    }

    return 0;
}
