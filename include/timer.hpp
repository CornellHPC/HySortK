#ifndef TIMER_H_
#define TIMER_H_

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <iomanip>

struct Timer
{
    bool isroot;
    MPI_Comm comm;
    double elapsed;

    Timer(MPI_Comm comm) : comm(comm), elapsed(0)
    {
        int myrank;
        MPI_Comm_rank(comm, &myrank);
        isroot = (myrank == 0);
    }

    void start()
    {
        MPI_Barrier(comm);
        elapsed = -MPI_Wtime();
    }

    void stop()
    {
        elapsed += MPI_Wtime();
    }

    void log(char const *label)
    {
        double maxtime, proctime;

        MPI_Reduce(&elapsed, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&elapsed, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (isroot)
        {
            std::cout << label << ":\n";
            std::cout << "    total time (user seconds): " << std::fixed << std::setprecision(3) << elapsed << "\n";
            std::cout << "    total cost (proc seconds): " << std::fixed << std::setprecision(3) << proctime << "\n" << std::endl;
            // std::cout << label << ": total time (user secs): " << std::fixed << std::setprecision(3) << elapsed << "\n"
                      // << label << ": total work (proc secs): " << std::fixed << std::setprecision(3) << proctime << "\n" << std::endl;
        }

        MPI_Barrier(comm);
    }

    void stop_and_log(char const *label)
    {
        stop();
        log(label);
    }
};

#endif