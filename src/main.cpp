#include<iostream>
#include<upcxx/upcxx.hpp>

#include "logger.hpp"
#include "timer.hpp"

int main(int argc, char **argv){
    upcxx::init();

    Logger log;

    Timer t;
    t.start();

    log() << "Hello from rank " << upcxx::rank_me() ;
    log.flush("Hello");

    /* some unbalanced work */
    int myrank = upcxx::rank_me();
    if(myrank == 2) {
        int temp = 0;
        for (int i = 0; i < 10000000; i++){
            temp = (temp + i) % 1763 + (i % 3) - i;
        }
        log()<<temp;
    }

    t.stop_and_log("Time is");

    upcxx::finalize();
    return 0;
}