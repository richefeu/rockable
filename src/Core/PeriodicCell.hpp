#ifndef PERIODICCELL_HPP
#define PERIODICCELL_HPP

#include "mat9.hpp"

class PeriodicCell {
 public:
  mat9r h;      ///< Matrix that hold the cell geometry (each column is a vector that defines a side)
  mat9r hinv;   ///< inverse of matrix h (recomputed each time h is updated)
  mat9r vh;     ///< Velocities of the collective DoF
  mat9r ah;     ///< Acceleration of the collective DoF
  mat9r Sig;    ///< Mean stress matrix over the cell
  double mass;  ///< Mass of the periodic cell (a somehow fictive data)

  PeriodicCell();
  PeriodicCell(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY, double ZZ);
  void Define(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY, double ZZ);
  void precomputeInverse();
  vec3r getBranchCorrection(const vec3r& ipos, const vec3r& jpos);
  void forceToStayInside(vec3r& s);
};

#endif /* end of include guard: PERIODICCELL_HPP */
