#ifndef PARADIS_SORT_HPP
#define PARADIS_SORT_HPP

/* This code is taken from github https://github.com/odanivan/simple_paradis/tree/master. it is working, maybe we should rewrite it. */

#include <omp.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

namespace paradis {

static const size_t RADIX_BITS = 8;
static const size_t MaxKisuu = 1 << RADIX_BITS;
static const size_t kRadixMask = (1 << RADIX_BITS) - 1;
static const size_t kRadixBin = 1 << RADIX_BITS;

namespace sort_utils {
template <typename T> inline void swap(T *a, T *b) {
  T t = *a;
  *a = *b;
  *b = t;
}

} // namespace sort_utils

template <typename T> long long determineDigitBucket(int stage, T num) {
  long long result = num.GetByte(stage);
  return result;
}

template <long long kth_byte, typename T>
inline void radix_sort_par(T *s, T *t, T *begin_itr, long long processes = 1) {
   // std::cout<<"kth_byte:"<<kth_byte<<std::endl;
  long long cnt[MaxKisuu] = {0};

  long long elenum = t - s;
  long long start = s - begin_itr;

  long long part = elenum / processes;
  long long res = elenum % processes;

  long long localHists[processes][MaxKisuu];
  long long gh[MaxKisuu], gt[MaxKisuu], starts[MaxKisuu];
  long long ph[processes][MaxKisuu];
  long long pt[processes][MaxKisuu];

  long long SumCi = elenum;
  long long pfp[processes + 1];
  long long var_p = processes;

#pragma omp parallel num_threads(processes)
  {
    long long th = omp_get_thread_num();
#pragma omp for
    for (long long i = 0; i < kRadixBin; ++i) {
      for (long long t = 0; t < processes; ++t)
        localHists[t][i] = 0;
    }
#pragma omp barrier
#pragma omp for
    for (long long i = start; i < start + elenum; i++) {
      long long digit = determineDigitBucket(kth_byte, *(begin_itr + i));
      localHists[th][digit]++;
    }
#pragma omp barrier
#pragma omp for
    for (long long i = 0; i < kRadixBin; i++) {
      for (long long j = 0; j < processes; j++) {
        cnt[i] += localHists[j][i];
      }
    }
#pragma omp barrier
#pragma omp single
    {
      gh[0] = start;
      gt[0] = gh[0] + cnt[0];
      starts[0] = gh[0];
    }

#pragma omp single
    for (long long i = 1; i < kRadixBin; i++) {
      gh[i] = gh[i - 1] + cnt[i - 1];
      gt[i] = gh[i] + cnt[i];
      starts[i] = gh[i];
    }
#pragma omp barrier

    while (SumCi != 0) {
#pragma omp for
      for (long long ii = 0; ii < processes; ii++) {
        long long pID = ii;
        for (long long i = 0; i < kRadixBin; i++) {
          part = (long long)(gt[i] - gh[i]) / (long long)var_p;
          res = (long long)(gt[i] - gh[i]) % (long long)(var_p);
          if (pID < var_p - 1) {
            ph[pID][i] = part * pID + gh[i];
            pt[pID][i] = part * (pID + 1LL) + gh[i];
          } else {
            ph[pID][i] = part * pID + gh[i];
            pt[pID][i] = part * (pID + 1LL) + gh[i] + res;
          }
        }

        for (long long i = 0; i < kRadixBin; i++) {
          long long head = ph[pID][i];
          while (head < pt[pID][i]) {
            T v = *(begin_itr + head);
            long long k = determineDigitBucket(kth_byte, v);
            while (k != i && ph[pID][k] < pt[pID][k]) {
              sort_utils::swap(&v, begin_itr + (long long)ph[pID][k]);
              ph[pID][k]++;
              k = determineDigitBucket(kth_byte, v);
            }
            if (k == i) {
              *(begin_itr + head) = *(begin_itr + ph[pID][i]);
              head++;
              *(begin_itr + ph[pID][i]) = v;
              ph[pID][i]++;
            } else {
              *(begin_itr + head) = v;
              head++;
            }
          }
        }
      }
#pragma omp single
      {
        SumCi = 0;
        long long pfpN = kRadixBin / var_p;
        long long pfpM = kRadixBin % var_p;
        pfp[0] = 0LL;
        long long pfpMR = 0LL;
        for (long long i = 1LL; i < var_p + 1LL; i++) {
          if (pfpMR < pfpM)
            pfpMR++;
          pfp[i] = i * pfpN + pfpMR;
        }
      }
#pragma omp barrier
#pragma omp for
      for (long long k = 0; k < processes; k++) {
        for (long long i = pfp[k]; i < pfp[k + 1]; i++) {
          long long tail = gt[i];
          {
            for (long long pID = 0; pID < processes; pID++) {
              long long head = ph[pID][i];
              while (head < pt[pID][i] && head < tail) {
                T v = *(begin_itr + head);
                head++;
                if (determineDigitBucket(kth_byte, v) != i) {
                  while (head <= tail) {
                    tail--;
                    T w = *(begin_itr + tail);
                    if (determineDigitBucket(kth_byte, w) == i) {
                      *(begin_itr + (head - 1)) = w;
                      *(begin_itr + tail) = v;
                      break;
                    }
                  }
                }
              }
            }
          }
          gh[i] = tail;
        }
      }
#pragma omp barrier
#pragma omp single
      {
        for (long long i = 0; i < kRadixBin; i++) {
          SumCi += (gt[i] - gh[i]);
        }
      }
#pragma omp barrier
    }
  }

  if (kth_byte > 0) {
#pragma omp parallel num_threads(processes)
#pragma omp single
    {
      for (long long i = 0; i < kRadixBin; i++) {
        long long nextStageThreads = 1;
        nextStageThreads =
            processes * (cnt[i] * (log(cnt[i]) / log(kRadixBin)) /
                         (elenum * (log(elenum) / log(kRadixBin))));
        if (cnt[i] > 64LL) {
#pragma omp task
          radix_sort_par<(kth_byte > 0 ? (kth_byte - 1) : 0)>(
              begin_itr + starts[i], begin_itr + (starts[i] + cnt[i]),
              begin_itr, std::max(nextStageThreads, 1LL));
        } else if (cnt[i] > 1) {
          std::sort(begin_itr + starts[i], begin_itr + (starts[i] + cnt[i]));
        }
      }
#pragma omp taskwait
    }
  }
}

template <typename T, size_t sz> void sort(T *begin, T *end, size_t thread_count) {
  
  radix_sort_par<sz - 1>(begin, end, begin, thread_count);

}

} // namespace paradis

#endif // PARADIS_SORT_HPP