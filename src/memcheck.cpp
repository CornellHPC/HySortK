/* Taken from https://hpcf.umbc.edu/general-productivity/checking-memory-usage/ */
#include "memcheck.hpp"

int get_memory_usage_kb(long* vmrss_kb, long* vmpeak_kb)
{
    /* Get the the current process' status file from the proc filesystem */
    FILE* procfile = fopen("/proc/self/status", "r");

    long to_read = 8192;
    char buffer[to_read];
    int read = fread(buffer, sizeof(char), to_read, procfile);
    fclose(procfile);

    short found_vmrss = 0;
    short found_vmpeak = 0;
    char* search_result;

    /* Look through proc status contents line by line */
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_vmrss == 0 || found_vmpeak == 0) )
    {
        search_result = strstr(line, "VmRSS:");
        if (search_result != NULL)
        {
            sscanf(line, "%*s %ld", vmrss_kb);
            found_vmrss = 1;
        }

        search_result = strstr(line, "VmHWM:");
        if (search_result != NULL)
        {
            sscanf(line, "%*s %ld", vmpeak_kb);
            found_vmpeak = 1;
        }

        line = strtok(NULL, delims);
    }

    return (found_vmrss == 1 && found_vmpeak == 1) ? 0 : 1;
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

int get_cluster_memory_usage_kb(long* vmrss_per_process, long* vmpeak_per_process, int root, int np)
{
    long vmrss_kb;
    long vmpeak_kb;
    int ret_code = get_memory_usage_kb(&vmrss_kb, &vmpeak_kb);

    if (ret_code != 0)
    {
        printf("Could not gather memory usage!\n");
        return ret_code;
    }

    MPI_Gather(&vmrss_kb, 1, MPI_UNSIGNED_LONG, 
        vmrss_per_process, 1, MPI_UNSIGNED_LONG, 
        root, MPI_COMM_WORLD);

    MPI_Gather(&vmpeak_kb, 1, MPI_UNSIGNED_LONG, 
        vmpeak_per_process, 1, MPI_UNSIGNED_LONG, 
        root, MPI_COMM_WORLD);

    return 0;
}


void print_mem_log(int nprocs, int myrank, std::string msg){
    long vmrss_per_process[nprocs];
    long vmpeak_per_process[nprocs];
    get_cluster_memory_usage_kb(vmrss_per_process, vmpeak_per_process, 0, nprocs);

    if (myrank == 0)
    {
        std::cout<<msg<<std::endl;
        
        uint64_t vmrss = 0;
        uint64_t vmpeak = 0;
        for (int i = 0; i < nprocs; i++)
        {
            vmrss += vmrss_per_process[i];
            vmpeak += vmpeak_per_process[i];
        }
        vmrss /= (1024 * 1024);
        vmpeak /= (1024 * 1024);

        printf("Total VmRSS = %6ld GB, VmHWM = %6ld GB\n", vmrss, vmpeak);
    }
}

void get_mem(int nprocs, int myrank, size_t& vmrss, size_t& vmsize){
    long vmrss_per_process[nprocs];
    long vmsize_per_process[nprocs];
    get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, nprocs);

    uint64_t vmrss_total = 0;
    uint64_t vmsize_total = 0;

    if (myrank == 0) {
        for (int i = 0; i < nprocs; i++)
        {
            vmrss_total += vmrss_per_process[i];
            vmsize_total += vmsize_per_process[i];
        }
        vmrss_total /= (1024 * 1024);
        vmsize_total /= (1024 * 1024);

        vmrss = vmrss_total;
        vmsize = vmsize_total;
    }
}