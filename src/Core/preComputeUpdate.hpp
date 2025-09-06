#pragma once

#include <cassert>
#include <omp.h>
#include <parallel/algorithm>
#include <vector>

// Reorganize activeInteractions by precomputing offset in Force for OpenMP.
class preComputeUpdate {
 public:
  // constructor
  preComputeUpdate() {}

  // Openmp version of a prefix sum
  inline void prefixsum_inplace(size_t* x, int N) {
    START_TIMER("preComputeUpdate::prefixsum");

    int nthreads;
#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    // size_t suma[nthreads + 1]; // vr: Variable Length Array (VLA), which is a feature borrowed from C99, but not part of standard C++
    std::vector<size_t> suma(nthreads + 1);
    suma[0] = 0;

#pragma omp parallel
    {
      const int ithread = omp_get_thread_num();
      size_t sum = 0;
#pragma omp for schedule(static)
      for (int i = 0; i < N; i++) {
        sum += x[i];
        x[i] = sum;
      }
      suma[ithread + 1] = sum;
#pragma omp barrier
      size_t offset = 0;
      for (int i = 0; i < (ithread + 1); i++) {
        offset += suma[i];
      }
#pragma omp for schedule(static)
      for (int i = 0; i < N; i++) {
        x[i] += offset;
      }
    }
  }

  inline int binarySearch(std::pair<size_t, size_t>* arr, int p, int r, size_t num) {
    if (p <= r) {
      int mid = (p + r) / 2;
      if (arr[mid].second == num) {
        return mid;
      }
      if (arr[mid].second > num) {
        return binarySearch(arr, p, mid - 1, num);
      }
      if (arr[mid].second < num) {
        return binarySearch(arr, mid + 1, r, num);
      }
    }
    return -1;
  }

  // resize vectors
  void resize(size_t nbOfInteractions, size_t nbOfParticles) {
    Ii.resize(nbOfInteractions);
    Ij.resize(nbOfInteractions);
    offset_i.resize(nbOfParticles);
    offset_j.resize(nbOfParticles);
    size_i.resize(nbOfParticles, 0);
    size_j.resize(nbOfParticles, 0);
  }

  void sort_Ij() {
    START_TIMER("preComputeUpdate::sort");
    __gnu_parallel::stable_sort(Ij.begin(), Ij.end(), [](std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) {
      return a.second < b.second;
    });
  }

  std::vector<std::pair<size_t, size_t>> Ii;  ///< Ii contains I->i for a given active interaction (first = id, second =
                                              ///< I->i). This vector is sorted in "accelerations" by construction.
  std::vector<std::pair<size_t, size_t>> Ij;  ///< Ij contains I->j for a given active interaction (first = id, second =
                                              ///< I->j). This vector is sorted in "accelerations".
  std::vector<size_t> offset_i;               ///< Store the precomputed offset on Ii for a given value of I->i
  std::vector<size_t> offset_j;               ///< Store the precomputed offset on Ij for a given value of I->j
  std::vector<size_t> size_i;  ///< Store the number of elements of the same value of Ii.second ie. size_i[id] = number
                               ///< of Elements in Ii.second equals to id.
  std::vector<size_t> size_j;  ///< Store the number of elements of the same value of Ij.second ie. size_j[id] = number
                               ///< of Elements in Ij.second equals to id.
};
