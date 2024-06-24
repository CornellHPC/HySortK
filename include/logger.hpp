#ifndef HYSORTK_LOGGER_H_
#define HYSORTK_LOGGER_H_

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>

namespace hysortk {

class Logger
{
    std::unique_ptr<std::ostringstream> logstream, rootstream;
    int myrank, nprocs;
    MPI_Comm comm;

    std::string prefix();

public:
    Logger(MPI_Comm comm);
    void flush(char const *label, int prank=-1);
    void flush(std::ostringstream& ss);
    void flush_root(char const *label);
    void flush(std::ostringstream& ss, int rank);
    std::ostringstream& operator()() { return *logstream; }
    std::string rankstr();
    std::string rankstr(int proc);
    static std::string readrangestr(size_t pos, size_t count);
};

} // namespace hysortk

#endif