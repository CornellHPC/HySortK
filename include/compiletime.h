#ifndef HYSORTK_COMPILE_TIME_H_
#define HYSORTK_COMPILE_TIME_H_

#include <limits>
#include <cstdint>

#ifndef KMER_SIZE
#error "KMER_SIZE must be defined"
#else
static_assert(2 < KMER_SIZE && KMER_SIZE < 96);
#ifdef SMER_SIZE
static_assert(0 < SMER_SIZE && SMER_SIZE <= KMER_SIZE);
#endif
#endif

#ifndef LOWER_KMER_FREQ
#error "LOWER_KMER_FREQ must be defined"
#elif !defined (UPPER_KMER_FREQ)
#error "UPPER_KMER_FREQ must be defined"
#else
static_assert(0 < LOWER_KMER_FREQ && LOWER_KMER_FREQ <= UPPER_KMER_FREQ && UPPER_KMER_FREQ <= std::numeric_limits<uint16_t>::max());
#endif

namespace hysortk {

typedef int32_t MPI_Count_t;
typedef int32_t MPI_Offset_t;
#define MPI_COUNT_TYPE MPI_INT

} // namespace hysortk

#endif
