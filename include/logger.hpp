#ifndef LOGGER_H_
#define LOGGER_H_

#include <upcxx/upcxx.hpp>
#include <iostream>

class Logger {
    std::unique_ptr<std::ostringstream> logstream, rootstream;
    upcxx::dist_object<std::string> str;
    int myrank, nprocs;

public:
    Logger();
    void flush(char const *label);
    void flush(std::ostringstream& ss);
    void flush(std::ostringstream& ss, int rank);
    std::ostringstream& operator()() { return *logstream; }
    static std::string readrangestr(size_t pos, size_t count);
};

#endif