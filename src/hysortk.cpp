#include<iostream>
#include<mpi.h>
#include<omp.h>
#include<fstream>

#include "logger.hpp"
#include "timer.hpp"
#include "fastaindex.hpp"
#include "dnabuffer.hpp"
#include "dnaseq.hpp"
#include "kmerops.hpp"
#include "compiletime.h"
#include "memcheck.hpp"
#include "hysortk.hpp"

namespace hysortk {

std::shared_ptr<DnaBuffer> read_dna_buffer(const std::string& fasta_fname, MPI_Comm comm){
    Timer timer(comm);
    timer.start();
    FastaIndex index(fasta_fname, comm);
    timer.stop_and_log("reading and scattering fasta file");

    timer.start();
    auto mydna = std::make_shared<DnaBuffer>(index.getmydna());
    timer.stop_and_log("reading and 2-bit encoding fasta sequences");

    return mydna;
}


std::unique_ptr<KmerListS> kmer_count(std::shared_ptr<DnaBuffer> mydna, MPI_Comm comm){
    /* start kmer counting */

    Timer timer(comm);
    Timer timer_kcount(MPI_COMM_WORLD);
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);


    // auto mydna = read_dna_buffer(fasta_fname, comm);

    size_t vm = 0;
    vm = get_mem_gb(nprocs, myrank, "VmRSS");
    if (myrank == 0){
        std::cout << "Memory usage (VmRSS): " << vm << " GB" << std::endl << std::endl;
    }

    timer_kcount.start();

    timer.start();
    auto tm = prepare_supermer(*mydna, comm);
    timer.stop_and_log("prepare_supermer");

    timer.start();
    exchange_supermer(tm, comm);
    timer.stop_and_log("exchange_supermer");

    timer.start();
    auto kmerlist = filter_kmer(tm, comm);
    timer.stop_and_log("filter_kmer");

    vm = get_mem_gb(nprocs, myrank, "VmHWM");
    if (myrank == 0){
        std::cout << "Peak memory usage (VmHWM): " << vm << " GB" << std::endl << std::endl;
    }

    timer_kcount.stop_and_log("Overall kmer counting (Excluding I/O)");

    return kmerlist;
}


void print_kmer_histogram(const KmerListS& kmerlist, MPI_Comm comm) {
    #if LOG_LEVEL >= 2

    Logger logger(comm);
    size_t maxcount = std::accumulate(kmerlist.cbegin(), kmerlist.cend(), 0, [](size_t cur, const auto& entry) { return std::max(cur, entry.cnt); });

    MPI_Allreduce(MPI_IN_PLACE, &maxcount, 1, MPI_INT, MPI_MAX, comm);

    std::vector<int> histo(maxcount+1, 0);

    for(size_t i = 0; i < kmerlist.size(); ++i)
    {
        int cnt = kmerlist[i].cnt;
        assert(cnt >= 1);
        histo[cnt]++;
    }

    MPI_Allreduce(MPI_IN_PLACE, histo.data(), maxcount+1, MPI_INT, MPI_SUM, comm);

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        std::cout << "#count\tnumkmers" << std::endl;

        for (int i = 1; i < histo.size(); ++i)
        {
            if (histo[i] > 0)
            {
                std::cout << i << "\t" << histo[i] << std::endl;
            }
        }
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
    #endif
}

void write_output_file(const KmerListS& kmerlist, const std::string& output_dir, MPI_Comm comm) {
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (myrank == 0)
    {
        std::cout << "Writing output files..." << std::endl;
    }

    std::string output_fname = output_dir + "/" + std::to_string(myrank) + ".out";


    std::ofstream ofs(output_fname);
    if (!ofs)
    {
        std::cerr << "Error: cannot open output file " << output_fname << std::endl;
        MPI_Abort(comm, 1);
    }

    for (size_t i = 0; i < kmerlist.size(); ++i)
    {
        ofs << kmerlist[i].kmer << "\t" << kmerlist[i].cnt << std::endl;
    }

}

} // namespace hysortk