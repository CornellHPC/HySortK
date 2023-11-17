#include <sys/sysinfo.h>
#include <iostream>

void mem_check(){
    struct sysinfo info;
    sysinfo(&info);
    long long total_ram = info.totalram * info.mem_unit;
    long long free_ram = info.freeram * info.mem_unit;
    long long available_ram = info.totalram * info.mem_unit;

    std::cout << "Total RAM: " << total_ram / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Free RAM: " << free_ram / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Available RAM: " << available_ram / (1024 * 1024) << " MB" << std::endl;
}