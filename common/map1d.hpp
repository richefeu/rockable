#ifndef MAP1D_HPP_1C0912F2
#define MAP1D_HPP_1C0912F2

template <typename T>
T map1d(T value, T istart, T istop, T ostart, T ostop) {
  return ostart + (ostop - ostart) * ((value - istart) / (istop - istart));
}

#endif /* end of include guard: MAP1D_HPP_1C0912F2 */
