#ifndef TIMER_H_
#define TIMER_H_

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <iomanip>

namespace hysortk {

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

    void local_start()
    {
        elapsed = -MPI_Wtime();
    }

    double local_stop()
    {
        elapsed += MPI_Wtime();
        return elapsed;
    }

    void log(char const *label)
    {
        double maxtime, proctime;

        MPI_Reduce(&elapsed, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&elapsed, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (isroot)
        {
            std::cout << label << ":\n";
            std::cout << "    total time (user seconds): " << std::fixed << std::setprecision(3) << maxtime << "\n";
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

    double get_elapsed() const
    {
        return elapsed + MPI_Wtime();
    }
};


struct TimerLocal
{
    double elapsed;

    TimerLocal() : elapsed(0)
    {
    }

    void start()
    {
        elapsed = -MPI_Wtime();
    }

    double stop()
    {
        elapsed += MPI_Wtime();
        return elapsed;
    }
};

} // namespace hysortk

#endif