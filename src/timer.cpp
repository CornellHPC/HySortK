#include "timer.hpp"

Timer::Timer()
{
    myrank = upcxx::rank_me();
    nprocs = upcxx::rank_n();
    start();
}

void Timer::start()
{
    st = std::chrono::high_resolution_clock::now();
}

void Timer::stop_and_log(char const *label)
{
    ed = std::chrono::high_resolution_clock::now();
    elapsed = ed - st;
    t = elapsed.count();
    upcxx::barrier();

    double tsum = upcxx::reduce_one(t, upcxx::op_fast_add, 0).wait();
    double tmax = upcxx::reduce_one(t, upcxx::op_fast_max, 0).wait();

    if (myrank == 0)
    {
            std::cout << "\n" << label << ":\n";
            std::cout << "    total time (user seconds): " << std::fixed << std::setprecision(3) << tmax << "\n";
            std::cout << "    total cost (proc seconds): " << std::fixed << std::setprecision(3) << tsum << "\n" << std::endl;
    }

    upcxx::barrier();

}