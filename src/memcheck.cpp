/* Taken from https://hpcf.umbc.edu/general-productivity/checking-memory-usage/ */
#include "memcheck.hpp"
#include <iostream>

int get_memory_usage_kb(long* vm_kb, std::string type)
{
    /* Get the the current process' status file from the proc filesystem */
    FILE* procfile = fopen("/proc/self/status", "r");

    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    short found_vm = 0;
    char* search_result;

    /* Look through proc status contents line by line */
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_vm == 0) )
    {
        search_result = strstr(line, (type+":").c_str());
        if (search_result != NULL)
        {
            sscanf(line, "%*s %ld", vm_kb);
            found_vm = 1;
        }
        line = strtok(NULL, delims);
    }

    return (found_vm == 1) ? 0 : 1;
}

int get_free_memory_kb(size_t* memfree_kb)
{
    /* Get the the current process' status file from the proc filesystem */
    FILE* procfile = fopen("/proc/meminfo", "r");

    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    short found_free = 0;
    char* search_result;

    /* Look through proc status contents line by line */
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_free==0) )
    {
        search_result = strstr(line, "MemFree:");
        if (search_result != NULL)
        {
            sscanf(line, "%*s %zu", memfree_kb);
            found_free = 1;
        }

        line = strtok(NULL, delims);
    }

    return (found_free) ? 0 : 1;
}

int get_cluster_memory_usage_kb(long* vm_per_process, int root, int np, std::string type)
{
    long vm_kb;
    int ret_code = get_memory_usage_kb(&vm_kb, type);

    if (ret_code != 0)
    {
        printf("Could not gather memory usage!\n");
        return ret_code;
    }

    MPI_Gather(&vm_kb, 1, MPI_UNSIGNED_LONG, 
        vm_per_process, 1, MPI_UNSIGNED_LONG, 
        root, MPI_COMM_WORLD);

    return 0;
}


size_t get_mem_gb(int nprocs, int myrank, std::string type){
    long vm_per_process[nprocs];
    if (get_cluster_memory_usage_kb(vm_per_process, 0, nprocs, type)!=0){
        std::cerr<<"Error in get_cluster_memory_usage_kb."<<std::endl;
        return 0;
    }

    size_t vm_total = 0;

    if (myrank == 0) {
        for (int i = 0; i < nprocs; i++)
        {
            vm_total += vm_per_process[i];
        }
        vm_total /= (1024 * 1024);
    }
    return vm_total;
}