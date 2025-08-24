#pragma once

#include "ProfilingTools.hpp"
#include <cassert>
#include <linkCells.hpp>
#include <mutex>

linkCells build_linkCells(AABB& Box, vec3r& CellMinSizes);

class Traversals {
  using Traversal = std::vector<AABB_Cell*>;
  using Wave = std::vector<Traversal>;

  bool exist = false;
  std::mutex* m_lock = nullptr;
  Traversal m_traversal;
  Wave m_waves;
  Wave m_block;

 public:
  Traversals() {}

  void setMutexSize(size_t size) {
    if (m_lock != nullptr) free(m_lock);
    m_lock = new std::mutex[size];
  }

  void buildTraversal(linkCells& gridCells) {
    START_TIMER("linkCells::buildTraversal");

    // regular method
    m_traversal.clear();
    auto& N = gridCells.N;
    auto& cells = gridCells.cells;

    for (size_t ix = 0; ix < N.x; ++ix) {
      for (size_t iy = 0; iy < N.y; ++iy) {
        for (size_t iz = 0; iz < N.z; ++iz) {
          if (cells[ix][iy][iz].bodies.size() > 0) {
            m_traversal.push_back(&cells[ix][iy][iz]);
          }
        }
      }
    }

    // wave method --> no block
    for (auto it : m_waves) it.clear();

    m_waves.clear();
    const size_t nwave = 27;
    m_waves.resize(nwave);  // do not forget the oversized bodies

    const vec3r wave[27] = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {0, 2, 0},
                            {1, 2, 0}, {2, 2, 0}, {0, 0, 1}, {1, 0, 1}, {2, 0, 1}, {0, 1, 1}, {1, 1, 1},
                            {2, 1, 1}, {0, 2, 1}, {1, 2, 1}, {2, 2, 1}, {0, 0, 2}, {1, 0, 2}, {2, 0, 2},
                            {0, 1, 2}, {1, 1, 2}, {2, 1, 2}, {0, 2, 2}, {1, 2, 2}, {2, 2, 2}};

#pragma omp parallel for
    for (size_t w = 0; w < nwave; w++) {
      auto source_point = wave[w];

      for (size_t ix = source_point.x; ix < N.x; ix += 3) {
        for (size_t iy = source_point.y; iy < N.y; iy += 3) {
          for (size_t iz = source_point.z; iz < N.z; iz += 3) {
            if (cells[ix][iy][iz].bodies.size() > 0) {
              m_waves[w].push_back(&cells[ix][iy][iz]);
            }
          }
        }
      }
    }

    // wave method --> with blocks
    for (auto it : m_block) it.clear();
    m_block.clear();
    constexpr size_t nblock = 8;
    constexpr size_t block_size = 2;
    assert(block_size >= 2);
    m_block.resize(nblock);

    const vec3r block_wave[nblock] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                      {0, 1, 1}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}};

#pragma omp parallel for
    for (size_t w = 0; w < nblock; w++) {
      auto source_point = block_wave[w];
      for (size_t ix = source_point.x * block_size; ix < N.x; ix += 2 * block_size) {
        for (size_t iy = source_point.y * block_size; iy < N.y; iy += 2 * block_size) {
          for (size_t iz = source_point.z * block_size; iz < N.z; iz += 2 * block_size) {
            // inner block
            for (size_t ix_block = 0; ix_block < block_size && (ix_block + ix) < N.x; ix_block++) {
              for (size_t iy_block = 0; iy_block < block_size && (iy_block + iy) < N.y; iy_block++) {
                for (size_t iz_block = 0; iz_block < block_size && (iz_block + iz) < N.z; iz_block++) {
                  if (cells[ix + ix_block][iy + iy_block][iz + iz_block].bodies.size() > 0) {
                    m_block[w].push_back(&cells[ix + ix_block][iy + iy_block][iz + iz_block]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Use it in an openemp parallel only.
  template <typename Func, typename I, typename... Args>
  inline void applyKernelTraversal(Traversal& a_traversal, Func& a_fun,
                                   std::vector<std::vector<I*>>& a_particles_interactions, Args&... a_args) {
#pragma omp for schedule(static)
    for (size_t idx = 0; idx < a_traversal.size(); idx++) {
      auto& particleIndex = a_traversal[idx]->bodies;
      for (auto& itParticles : particleIndex) {
        auto& interactions = a_particles_interactions[itParticles];
        for (auto& interaction : interactions) {
          a_fun(*interaction, a_args...);
        }
      }
    }
  }

  // Oversized_bodies are not procceed by this function
  template <typename Func, typename I, typename... Args>
  inline void applyKernelMutexes(Func& a_fun, std::vector<std::vector<I*>>& a_particles_interactions, Args&... a_args) {
    START_TIMER("linkCells::apply_kernel_mutexes");
#pragma omp parallel
    applyKernelTraversal(m_traversal, a_fun, a_particles_interactions, m_lock, a_args...);
  }

  // oversized_bodies are not procceed by this function
  template <typename Func, typename I, typename... Args>
  inline void applyKernel(Func& a_fun, std::vector<std::vector<I*>>& a_particles_interactions, Args&... a_args) {
    START_TIMER("linkCells::apply_kernel_wave");
#pragma omp parallel
    for (size_t w = 0; w < m_waves.size(); w++) {
      applyKernelTraversal(m_waves[w], a_fun, a_particles_interactions, a_args...);
    }
  }

  // oversized_bodies are not procceed by this function
  template <typename Func, typename I, typename... Args>
  inline void applyKernelBlock(Func& a_fun, std::vector<std::vector<I*>>& a_particles_interactions, Args&... a_args) {
    START_TIMER("linkCells::apply_kernel_block");
#pragma omp parallel
    for (size_t w = 0; w < m_block.size(); w++) {
      applyKernelTraversal(m_block[w], a_fun, a_particles_interactions, a_args...);
    }
  }

  // accessor
  inline std::vector<size_t>& get_oversized_bodies_indexes(linkCells& lc) {
    auto& ret = lc.oversized_bodies.bodies;
    return ret;
  }

  // accessor
  inline std::mutex* get_mutexes() {
    auto ret = m_lock;
    return ret;
  }

  // proceed oversized particles (with a huge number of interactions per particle)
  template <typename Func, typename I, typename... Args>
  inline void apply_kernel_oversized_bodies_only(linkCells& lc, Func& a_fun,
                                                 std::vector<std::vector<I*>>& a_particles_interactions,
                                                 Args&... a_args) {
    START_TIMER("linkCells::apply_kernel_oversized_bodies_only");
    auto& indexes = get_oversized_bodies_indexes(lc);
    std::mutex* mutexes = get_mutexes();
    using field = std::vector<I*>;

#pragma omp parallel
    for (auto it : indexes) {
      const size_t size = a_particles_interactions[it].size();
      field& vecInteractions = a_particles_interactions[it];
#pragma omp for nowait
      for (size_t i = 0; i < size; i++) {
        I& Interaction = *vecInteractions[i];
        a_fun(Interaction, mutexes, a_args...);
      }
    }
  }
};
