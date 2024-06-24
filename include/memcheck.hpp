# ifndef HYSORTK_MEMCHECK_HPP
# define HYSORTK_MEMCHECK_HPP

#include <mpi.h>
#include <memory.h>
#include <cstdio>
#include <iostream>
#include <string>

namespace hysortk {

int get_free_memory_kb(size_t* memfree_kb);
size_t get_mem_gb(int nprocs, int myrank, std::string type);

} // namespace hysortk

# endif