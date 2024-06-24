#ifndef HYSORTK_DNABUFFER_H_
#define HYSORTK_DNABUFFER_H_

#include "dnaseq.hpp"
#include <memory>
#include <numeric>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <string>

namespace hysortk {

class DnaBuffer
{
public:
    DnaBuffer(size_t bufsize) : bufhead(0), bufsize(bufsize), buf(new uint8_t[bufsize]) {}
    DnaBuffer(size_t bufsize, size_t numreads, uint8_t *buf, const size_t *readlens);
    DnaBuffer(const DnaBuffer& other) : bufhead(other.bufhead), bufsize(other.bufsize), buf(new uint8_t[bufsize]) {
        std::copy(other.buf, other.buf + bufsize, buf);
        sequences.clear();
        size_t offset = 0;
        for (size_t i = 0; i < other.size(); ++i) {
            sequences.emplace_back(other[i].size(), buf + offset);
            offset += other[i].numbytes();
        }
    }

    void push_back(char const *s, size_t len);
    size_t size() const { return sequences.size(); }
    size_t getbufsize() const { return bufsize; }
    size_t getrangebufsize(size_t start, size_t count) const;
    const uint8_t* getbufoffset(size_t i) const { return sequences[i].data(); }
    const DnaSeq& operator[](size_t i) const { return sequences[i]; }

    std::string getasciifilecontents() const;

    static size_t computebufsize(const std::vector<size_t>& seqlens);

    ~DnaBuffer() { delete[] buf; }

private:
    size_t bufhead;
    const size_t bufsize;
    uint8_t *buf;
    std::vector<DnaSeq> sequences;
};

} // namespace hysortk

#endif