#include<iostream>
#include<upcxx/upcxx.hpp>

#include "logger.hpp"
#include "timer.hpp"
#include "fastaindex.hpp"
#include "dnabuffer.hpp"
#include "dnaseq.hpp"

std::string fasta_fname;


int main(int argc, char **argv){
    upcxx::init();

    Logger log;
    Timer timer;
    std::ostringstream ss;

    if (argc < 2){
        std::cerr << "Usage: " << argv[0] << " <fasta file>" << std::endl;
        exit(1);
    }

    fasta_fname = argv[1];

    timer.start();
    FastaIndex index(fasta_fname);
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

    upcxx::finalize();
    return 0;
}