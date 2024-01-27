# ifndef MEMCHECK_HPP
# define MEMCHECK_HPP

#include <mpi.h>
#include <memory.h>
#include <cstdio>
#include <iostream>
#include <string>

int get_cluster_memory_usage_kb(long* vmrss_per_process, long* vmpeak_per_process, int root, int np);
int get_global_memory_usage_kb(long* global_vmrss, long* global_vmpeak, int np);
void print_mem_log(int nprocs, int myrank, std::string msg=std::string());

# endif