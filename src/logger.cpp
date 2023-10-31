#include "logger.hpp"
#include <iostream>
#include <numeric>
#include <upcxx/upcxx.hpp>

Logger::Logger() : logstream(new std::ostringstream())
{
    myrank = upcxx::rank_me();
    nprocs = upcxx::rank_n();
    str = upcxx::dist_object<std::string>("");

}


void Logger::flush(std::ostringstream& ss)
{
    flush(ss.str().c_str());
    ss.clear();
    ss.str("");
}

void Logger::flush(std::ostringstream& ss, int rank)
{
    if (rank == myrank)
    {
        std::cout << ss.str().c_str() << std::endl;
    }
    ss.clear();
    ss.str("");

    upcxx::barrier();
}

void Logger::flush(char const *label)
{
    int myrank = upcxx::rank_me();
    int nprocs = upcxx::rank_n();

    str = logstream->str();
    logstream.reset(new std::ostringstream());

    upcxx::barrier();

    if (myrank == 0)
    {
        std::string slabel(label);
        std::string banner;
        banner.assign(slabel.size(), '=');
        std::cout << slabel << "\n" << banner << "\n" << std::endl;

        for (int i = 0; i < nprocs; i++)
        {
            std::string s = str.fetch(i).wait();
            std::cout << "rank[" << i+1 << "/" << nprocs << "] :: " << s << std::endl;
        }
    }

    upcxx::barrier();

}

std::string Logger::readrangestr(size_t pos, size_t count)
{
    std::ostringstream ss;
    ss << "[" << pos << ".." << (pos+count-1) << "] (" << count << " reads)";
    return ss.str();
}
