#ifndef MMER_H_
#define MMER_H_

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>
#include "compiletime.h"
#include "dnaseq.hpp"
#include "hashfuncs.hpp"
#include "dnaseq.hpp"

#define MAX_SUPERMER_LEN 250

template <int MLONGS>
class Mmer
{
public:

    static_assert(MLONGS != 0);

    static constexpr int NBYTES = 8 * MLONGS;

    typedef std::array<uint64_t, MLONGS> MERARR;
    typedef std::array<uint8_t,  NBYTES> BYTEARR;

    Mmer();
    Mmer(const DnaSeq& s);
    Mmer(char const *s);
    Mmer(const void *mem);
    Mmer(const Mmer& o);

    Mmer& operator=(Mmer o);

    std::string GetString() const;

    bool operator<(const Mmer& o) const;
    bool operator==(const Mmer& o) const;
    bool operator!=(const Mmer& o) const;

    Mmer GetExtension(int code) const;
    Mmer GetTwin() const;
    Mmer GetRep() const;

    uint64_t GetHash() const;
    const void* GetBytes() const { return reinterpret_cast<const void*>(longs.data()); }
    int getByte(int &i) const {
        return bytes[i];
    }

    void CopyDataInto(void *mem) const { std::memcpy(mem, longs.data(), NBYTES); }
    void CopyDataFrom(const void *mem) { std::memcpy(longs.data(), mem, NBYTES); }

    static std::vector<Mmer> GetMmers(const DnaSeq& s);
    static std::vector<Mmer> GetRepMmers(const DnaSeq& s);

    template <int N>
    friend std::ostream& operator<<(std::ostream& os, const Mmer<N>& kmer);

private:

    union { MERARR  longs;
            BYTEARR bytes; };

    void set_kmer(const DnaSeq& s);
    void set_kmer(char const *s, bool const revcomp = false);
};

template <int MLONGS>
std::ostream& operator<<(std::ostream& os, const Mmer<MLONGS>& kmer)
{
    os << MINIMIZER_SIZE << "-mer(" << kmer.GetString() << ")";
    return os;
}

namespace std
{
    template <int MLONGS> struct hash<Mmer<MLONGS>>
    {
        size_t operator()(const Mmer<MLONGS>& kmer) const
        {
            auto myhash = kmer.GetHash();
            return myhash;
        }
    };

    template <int MLONGS> struct less<Mmer<MLONGS>>
    {
        bool operator()(const Mmer<MLONGS>& k1, const Mmer<MLONGS>& k2) const
        {
            return k1 < k2;
        }
    };
}

/*
static uint64_t tetramer_twin(const uint8_t code)
{
    static const uint8_t tetramer_lookup_code[256] =
    {
        0xff,0xbf,0x7f,0x3f,0xef,0xaf,0x6f,0x2f,0xdf,0x9f,0x5f,0x1f,0xcf,0x8f,0x4f,0xf,
        0xfb,0xbb,0x7b,0x3b,0xeb,0xab,0x6b,0x2b,0xdb,0x9b,0x5b,0x1b,0xcb,0x8b,0x4b,0xb,
        0xf7,0xb7,0x77,0x37,0xe7,0xa7,0x67,0x27,0xd7,0x97,0x57,0x17,0xc7,0x87,0x47,0x7,
        0xf3,0xb3,0x73,0x33,0xe3,0xa3,0x63,0x23,0xd3,0x93,0x53,0x13,0xc3,0x83,0x43,0x3,
        0xfe,0xbe,0x7e,0x3e,0xee,0xae,0x6e,0x2e,0xde,0x9e,0x5e,0x1e,0xce,0x8e,0x4e,0xe,
        0xfa,0xba,0x7a,0x3a,0xea,0xaa,0x6a,0x2a,0xda,0x9a,0x5a,0x1a,0xca,0x8a,0x4a,0xa,
        0xf6,0xb6,0x76,0x36,0xe6,0xa6,0x66,0x26,0xd6,0x96,0x56,0x16,0xc6,0x86,0x46,0x6,
        0xf2,0xb2,0x72,0x32,0xe2,0xa2,0x62,0x22,0xd2,0x92,0x52,0x12,0xc2,0x82,0x42,0x2,
        0xfd,0xbd,0x7d,0x3d,0xed,0xad,0x6d,0x2d,0xdd,0x9d,0x5d,0x1d,0xcd,0x8d,0x4d,0xd,
        0xf9,0xb9,0x79,0x39,0xe9,0xa9,0x69,0x29,0xd9,0x99,0x59,0x19,0xc9,0x89,0x49,0x9,
        0xf5,0xb5,0x75,0x35,0xe5,0xa5,0x65,0x25,0xd5,0x95,0x55,0x15,0xc5,0x85,0x45,0x5,
        0xf1,0xb1,0x71,0x31,0xe1,0xa1,0x61,0x21,0xd1,0x91,0x51,0x11,0xc1,0x81,0x41,0x1,
        0xfc,0xbc,0x7c,0x3c,0xec,0xac,0x6c,0x2c,0xdc,0x9c,0x5c,0x1c,0xcc,0x8c,0x4c,0xc,
        0xf8,0xb8,0x78,0x38,0xe8,0xa8,0x68,0x28,0xd8,0x98,0x58,0x18,0xc8,0x88,0x48,0x8,
        0xf4,0xb4,0x74,0x34,0xe4,0xa4,0x64,0x24,0xd4,0x94,0x54,0x14,0xc4,0x84,0x44,0x4,
        0xf0,0xb0,0x70,0x30,0xe0,0xa0,0x60,0x20,0xd0,0x90,0x50,0x10,0xc0,0x80,0x40,0x0
    };

    return static_cast<uint64_t>(tetramer_lookup_code[code]);
}
*/

template <int MLONGS>
Mmer<MLONGS>::Mmer() : longs{} {}

template <int MLONGS>
Mmer<MLONGS>::Mmer(const DnaSeq& s) : Mmer() { set_kmer(s); }

template <int MLONGS>
Mmer<MLONGS>::Mmer(char const *s) : Mmer() { set_kmer(s); }

template <int MLONGS>
Mmer<MLONGS>::Mmer(const Mmer& o) : longs(o.longs) {}

template <int MLONGS>
Mmer<MLONGS>::Mmer(const void *mem) : Mmer() { CopyDataFrom(mem); }

template <int MLONGS>
std::string Mmer<MLONGS>::GetString() const
{
    std::string s(MINIMIZER_SIZE, '\0');

    int i, j, l;

    for (i = 0; i < MINIMIZER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        s[i] = "ACGT"[(longs[l] >> (2 * (31 - j)))&3];
    }

    return s;
}

template <int MLONGS>
void Mmer<MLONGS>::set_kmer(const DnaSeq& s)
{
    int i, j, l, idx;
    uint64_t code;

    /*
     * Warning: set_kmer assumes that longs/bytes have already
     * been completely zeroed out.
     */

    for (i = 0; i < MINIMIZER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        code = static_cast<uint64_t>(s[i]);

        longs[l] |= (code << (2 * (31 - j)));
    }
}

template <int MLONGS>
void Mmer<MLONGS>::set_kmer(char const *s, bool const revcomp)
{
    int i, j, l, idx;
    uint64_t code;

    /*
     * Warning: set_kmer assumes that longs/bytes have already
     * been completely zeroed out.
     */

    for (i = 0; i < MINIMIZER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        idx = revcomp? MINIMIZER_SIZE - i - 1 : i;
        code = static_cast<uint64_t>(DnaSeq::getcharcode(s[idx]));

        longs[l] |= ((revcomp? 3 - code : code) << (2 * (31 - j)));
    }
}
template <int MLONGS>
Mmer<MLONGS>& Mmer<MLONGS>::operator=(Mmer o)
{
    std::swap(longs, o.longs);
    return *this;
}

template <int MLONGS>
bool Mmer<MLONGS>::operator<(const Mmer& o) const
{
    for (int i = 0; i < MLONGS; ++i)
    {
        if (longs[i] < o.longs[i])
            return true;

        if (longs[i] > o.longs[i])
            return false;
    }

    return false;
}

template <int MLONGS>
bool Mmer<MLONGS>::operator==(const Mmer& o) const
{
    for (int i = 0; i < MLONGS; ++i)
        if (longs[i] != o.longs[i])
            return false;

    return true;
}

template <int MLONGS>
bool Mmer<MLONGS>::operator!=(const Mmer& o) const
{
    return !(*this == o);
}

template <int MLONGS>
Mmer<MLONGS> Mmer<MLONGS>::GetExtension(int code) const
{
    Mmer ext;

    ext.longs[0] = longs[0] << 2;

    for (uint64_t i = 1; i < MLONGS; ++i)
    {
        ext.longs[i-1] |= ((longs[i] >> 62) & 0x3);
        ext.longs[i] = longs[i] << 2;
    }

    ext.longs[MLONGS-1] |= (static_cast<uint64_t>(code) << (2 * (32 - (MINIMIZER_SIZE%32))));

    return ext;
}

template <int MLONGS>
Mmer<MLONGS> Mmer<MLONGS>::GetTwin() const
{
    Mmer twin;

    /* unroll */
    for (int l = 0; l < MLONGS; ++l)
    {
        uint64_t longmer = longs[l];

        /* unroll */
        for (uint64_t i = 0; i < 64; i += 8)
        {
            uint8_t bytemer = (longmer >> i) & 0xff;
            uint64_t revcomp_bytemer = tetramer_twin(bytemer);
            twin.longs[MLONGS-1-l] |= (revcomp_bytemer << (56 - i));
        }
    }

    uint64_t shift = MINIMIZER_SIZE % 32? 2 * (32 - (MINIMIZER_SIZE % 32)) : 0ULL;
    uint64_t mask = MINIMIZER_SIZE % 32? ((1ULL << shift) - 1) << (64 - shift) : 0ULL;

    twin.longs[0] <<= shift;

    for (uint64_t i = 1; i < MLONGS; ++i)
    {
        twin.longs[i-1] |= (twin.longs[i] & mask) >> (64 - shift);
        twin.longs[i] <<= shift;
    }

    return twin;
}

template <int MLONGS>
Mmer<MLONGS> Mmer<MLONGS>::GetRep() const
{
    Mmer twin = GetTwin();
    return twin < *this? twin : *this;
}

template <int MLONGS>
uint64_t Mmer<MLONGS>::GetHash() const
{
    uint64_t h;
    murmurhash3_64(longs.data(), NBYTES, &h);
    return h;
}

template <int MLONGS>
std::vector<Mmer<MLONGS>> Mmer<MLONGS>::GetMmers(const DnaSeq& s)
{
    int l = s.size();
    int num_kmers = l - MINIMIZER_SIZE + 1;

    if (num_kmers <= 0) return std::vector<Mmer>();

    std::vector<Mmer> kmers;

    kmers.reserve(num_kmers);
    kmers.emplace_back(s);

    for (int i = 1; i < num_kmers; ++i)
    {
        kmers.push_back(kmers.back().GetExtension(s[i+MINIMIZER_SIZE-1]));
    }

    return kmers;
}

template <int MLONGS>
std::vector<Mmer<MLONGS>> Mmer<MLONGS>::GetRepMmers(const DnaSeq& s)
{
    auto kmers = GetMmers(s);
    std::transform(kmers.begin(), kmers.end(), kmers.begin(), [](const Mmer& kmer) { return kmer.GetRep(); });
    return kmers;
}


using TMmer = typename std::conditional<(MINIMIZER_SIZE <= 32), Mmer<1>,
              typename std::conditional<(MINIMIZER_SIZE <= 64), Mmer<2>,
              typename std::conditional<(MINIMIZER_SIZE <= 96), Mmer<3>, Mmer<0>>::type>::type>::type;

#endif
