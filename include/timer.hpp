#ifndef TIMER_H_
#define TIMER_H_

#include <upcxx/upcxx.hpp>
#include <iostream>
#include <chrono>

class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> st, ed;
    std::chrono::duration<double> elapsed;
    double t;
    int myrank, nprocs;

public:
    Timer();
    void start();
    void stop_and_log(char const *label);
};

#endif