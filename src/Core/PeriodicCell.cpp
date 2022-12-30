#include "PeriodicCell.hpp"

PeriodicCell::PeriodicCell() : h(), vh(), ah() {}

/**
 * @brief Construct a new Periodic Cell:: Periodic Cell object. Defines each component of matrix h (cell geometry) 
 */
PeriodicCell::PeriodicCell(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY,
                           double ZZ)
    : h(XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ), vh(), ah() {
}

void PeriodicCell::precomputeInverse() {
  hinv = h.get_inverse();
}

/**
 * @brief Get the branch correction to account for the periodicity.
 * 
 * @param ipos real coordinates of the first body (i)
 * @param jpos real coordinates of the second body (j)
 * @return vec3r a vector (real coordinates) to be add to the position of the second body. If this vector is non-zero,
 *               then adding it to the position of j will be its image
 */
vec3r PeriodicCell::getBranchCorrection(const vec3r& ipos, const vec3r& jpos) {
  vec3r rij = jpos - ipos;
  vec3r sij = hinv * rij;
  sij.x -= floor(sij.x + 0.5);
  sij.y -= floor(sij.y + 0.5);
  sij.z -= floor(sij.z + 0.5);
  return (h * sij - rij);
}

/**
 * @brief restricts the reduced coordinates to be in range [0 1]
 * 
 * @param s the reduced coordinate to be modified 
 */
void PeriodicCell::forceToStayInside(vec3r& s) {
  s.x -= floor(s.x);
  s.y -= floor(s.y);
  s.z -= floor(s.z);
}
