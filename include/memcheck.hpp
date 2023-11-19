# ifndef MEMCHECK_HPP
# define MEMCHECK_HPP

#include <mpi.h>
#include <memory.h>
#include <cstdio>
#include <iostream>
#include <string>

int get_cluster_memory_usage_kb(long* vmrss_per_process, long* vmsize_per_process, int root, int np);
int get_global_memory_usage_kb(long* global_vmrss, long* global_vmsize, int np);
void print_mem_log(int nprocs, int myrank, std::string msg=std::string());
void get_mem(int nprocs, int myrank, size_t& vmrss, size_t& vmsize);

class memchecker{

    size_t vmrss;
    size_t vmsize;
    int nprocs;
    int myrank;

public:
    void init() {
        get_mem(nprocs, myrank, vmrss, vmsize);
    }

    memchecker(int nprocs, int myrank) : nprocs(nprocs), myrank(myrank) {
        init();
    }

    void log(std::string msg=std::string()) {
        size_t vmrss_new;
        size_t vmsize_new;
        get_mem(nprocs, myrank, vmrss_new, vmsize_new);

        if (myrank == 0) {
            std::cout << msg << ": " << std::endl;
            std::cout << "VmRSS: " << vmrss_new - vmrss << " GB" << std::endl;
            std::cout << "VmSize: " << vmsize_new - vmsize << " GB" << std::endl;
            std::cout << std::endl;
        }

    }
};

# endif