#ifndef PERIODIZOR_HPP
#define PERIODIZOR_HPP

#include <cmath>

#include "mat9.hpp"

class Periodizor {

  mat9r cell;
  mat9r cellUnitDir;
  vec3r invLength;
  vec3r length;

public:
  Periodizor() {
    vec3r one(1.0, 1.0, 1.0);
    update(one);
  }
  Periodizor(vec3r &cellAxisAligned) { update(cellAxisAligned); }
  Periodizor(mat9r &cell) { update(cell); }

  void update(vec3r &cellAxisAligned) {
    cell.xx = cellAxisAligned.x;
    cell.xy = 0.0;
    cell.xz = 0.0;
    cell.yx = 0.0;
    cell.yy = cellAxisAligned.y;
    cell.yz = 0.0;
    cell.zx = 0.0;
    cell.zy = 0.0;
    cell.zz = cellAxisAligned.z;
    length = cellAxisAligned;
    invLength.x = 1.0 / length.x;
    invLength.y = 1.0 / length.y;
    invLength.z = 1.0 / length.z;
    cellUnitDir.xx = 1.0;
    cellUnitDir.xy = 0.0;
    cellUnitDir.xz = 0.0;
    cellUnitDir.yx = 0.0;
    cellUnitDir.yy = 1.0;
    cellUnitDir.yz = 0.0;
    cellUnitDir.zx = 0.0;
    cellUnitDir.zy = 0.0;
    cellUnitDir.zz = 1.0;
  }

  void update(mat9r &cell_) {
    cell = cell_;
    length.x = sqrt(cell.xx * cell.xx + cell.yx * cell.yx + cell.zx * cell.zx);
    length.y = sqrt(cell.xy * cell.xy + cell.yy * cell.yy + cell.zy * cell.zy);
    length.z = sqrt(cell.xz * cell.xz + cell.yz * cell.yz + cell.zz * cell.zz);
    invLength.x = 1.0 / length.x;
    invLength.y = 1.0 / length.y;
    invLength.z = 1.0 / length.z;
    cellUnitDir.xx = cell.xx * invLength.x;
    cellUnitDir.yx = cell.yx * invLength.x;
    cellUnitDir.zx = cell.zx * invLength.x;
    cellUnitDir.xy = cell.xy * invLength.y;
    cellUnitDir.yy = cell.yy * invLength.y;
    cellUnitDir.zy = cell.zy * invLength.y;
    cellUnitDir.xz = cell.xz * invLength.z;
    cellUnitDir.yz = cell.yz * invLength.z;
    cellUnitDir.zz = cell.zz * invLength.z;
  }

  void tune(vec3r &vecToBeTuned) {
    vec3r shift(
        (cellUnitDir.xx * vecToBeTuned.x + cellUnitDir.yx * vecToBeTuned.y +
         cellUnitDir.zx * vecToBeTuned.z) *
            invLength.x,
        (cellUnitDir.xy * vecToBeTuned.x + cellUnitDir.yy * vecToBeTuned.y +
         cellUnitDir.zy * vecToBeTuned.z) *
            invLength.y,
        (cellUnitDir.xz * vecToBeTuned.x + cellUnitDir.yz * vecToBeTuned.y +
         cellUnitDir.zz * vecToBeTuned.z) *
            invLength.z);

    shift.x = std::round(shift.x);
    shift.y = std::round(shift.y);
    shift.z = std::round(shift.z);

    vec3r vshift;
    if (shift.x != 0.0)
      vshift += cell.get_xcol() * shift.x;
    if (shift.y != 0.0)
      vshift += cell.get_ycol() * shift.y;
    if (shift.z != 0.0)
      vshift += cell.get_zcol() * shift.z;

    vecToBeTuned -= vshift;
  }

  void tuneAxisAligned(vec3r &vecToBeTuned) {
    if (vecToBeTuned.x != 0.0) {
      double shift = vecToBeTuned.x * invLength.x;
      shift = std::round(shift) * length.x;
      vecToBeTuned.x -= shift;
    }
    if (vecToBeTuned.y != 0.0) {
      double shift = vecToBeTuned.y * invLength.y;
      shift = std::round(shift) * length.y;
      vecToBeTuned.y -= shift;
    }
    if (vecToBeTuned.z != 0.0) {
      double shift = vecToBeTuned.z * invLength.z;
      shift = std::round(shift) * length.z;
      vecToBeTuned.z -= shift;
    }
  }
};

#if 1
int main() {

  mat9r h(10, 0, 0, 5, 10, 0, 0, 0, 1);

  Periodizor P(h);

  // vec3r pp(10., 10., 10.);
  // Periodizor P(pp);

  vec3r p(8., 13., 0.0);
  P.tune(p);
  std::cout << p << '\n';

  return 0;
}
#endif

#endif /* end of include guard: PERIODIZOR_HPP */
