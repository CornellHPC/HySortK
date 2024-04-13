#include "logger.hpp"
#include <iostream>
#include <numeric>

Logger::Logger(MPI_Comm comm) : logstream(new std::ostringstream()), comm(comm)
{
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
}

std::string Logger::readrangestr(size_t pos, size_t count)
{
    std::ostringstream ss;
    ss << "[" << pos << ".." << (pos+count-1) << "] (" << count << " reads)";
    return ss.str();
}

std::string Logger::rankstr(int proc)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);

    std::ostringstream ss;
    ss << "rank[" << proc+1 << "/" << nprocs << "]";
    return ss.str();
}

std::string Logger::rankstr()
{
    static bool initialized = false;
    static std::string name;

    if (!initialized)
    {
        int myrank;
        int nprocs;
        MPI_Comm_rank(comm, &myrank);
        MPI_Comm_size(comm, &nprocs);

        std::ostringstream ss;
        ss << "rank[" << myrank+1 << "/" << nprocs << "]";
        name = ss.str();
        initialized = true;
    }

    return name;
}

std::string Logger::prefix()
{
    static bool initialized = false;
    static std::string name;

    if (!initialized)
    {
        int myrank;
        int nprocs;
        MPI_Comm_rank(comm, &myrank);
        MPI_Comm_size(comm, &nprocs);

        std::ostringstream ss;
        ss << "rank[" << myrank+1 << "/" << nprocs << "] :: ";
        name = ss.str();
        initialized = true;
    }

    return name;
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
        std::cout << ss.str() << std::endl;
    }

    ss.clear();
    ss.str("");

    MPI_Barrier(comm);
}

void Logger::flush_root(char const *label)
{

    std::string mylog = logstream->str();

    if (myrank == 0)
    {
        std::string slabel(label);
        std::cout<< slabel << std::endl;
        std::cout << mylog << std::endl << std::endl;
    }

    logstream.reset(new std::ostringstream());

    MPI_Barrier(comm);
}

void Logger::flush(char const *label, int prank)
{
    int myrank;
    int nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<int> recvcnt, displs;
    std::vector<char> recvbuf;

    std::string mylog = logstream->str();
    logstream.reset(new std::ostringstream());

    int sendcnt = mylog.size();
    int totrecv;

    if (!myrank) recvcnt.resize(nprocs);

    MPI_Gather(&sendcnt, 1, MPI_INT, recvcnt.data(), 1, MPI_INT, 0, comm);

    if (!myrank)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        std::partial_sum(recvcnt.begin(), recvcnt.end()-1, displs.begin()+1);
        recvbuf.resize(recvcnt.back() + displs.back());
    }

    MPI_Gatherv(mylog.c_str(), sendcnt, MPI_CHAR, recvbuf.data(), recvcnt.data(), displs.data(), MPI_CHAR, 0, comm);

    if (!myrank)
    {
        std::string slabel(label);
        std::string banner;
        banner.assign(slabel.size(), '=');
        std::cout << slabel << "\n" << banner << "\n" << std::endl;

        char const *buf = recvbuf.data();

        if (prank == -1) {
            for (int i = 0; i < nprocs; ++i)
            {
                std::cout << "rank[" << i+1 << "/" << nprocs << "] :: " << std::string(buf + displs[i], recvcnt[i]) << "\n";
            }
            std::cout << std::endl;
        } else {
            std::cout << "rank[-/-] :: " << std::string(buf + displs[0], recvcnt[0]) << "\n\n";
        }
    }

    MPI_Barrier(comm);
}


