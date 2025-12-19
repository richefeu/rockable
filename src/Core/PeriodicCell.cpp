//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "PeriodicCell.hpp"

PeriodicCell::PeriodicCell() : h(), vh(), ah() {}

/// @brief Construct a new Periodic Cell:: Periodic Cell object. Defines each component of matrix h (cell geometry)
///
PeriodicCell::PeriodicCell(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY,
                           double ZZ)
    : h(XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ), vh(), ah() {}

void PeriodicCell::precomputeInverse() { hinv = h.get_inverse(); }

/// @brief Get the branch correction to account for the periodicity.
///
/// @param ipos real coordinates of the first body (i)
/// @param jpos real coordinates of the second body (j)
/// @return vec3r a vector (real coordinates) to be add to the position of the second body. If this vector is non-zero,
///               then adding it to the position of j will be its image
///
vec3r PeriodicCell::getBranchCorrection(const vec3r& ipos, const vec3r& jpos) {
  vec3r rij = ipos - jpos;  // the branch vector from i to j
  vec3r sij = hinv * rij;   // same vector in reduced coordinates

  sij.x = std::round(sij.x);
  sij.y = std::round(sij.y);
  sij.z = std::round(sij.z);

  return (h * sij);  // back to real coordinates
}

/// @brief restricts the reduced coordinates to be in range [0 1]
///
/// @param s the reduced coordinate to be modified
///
void PeriodicCell::forceToStayInside(vec3r& s) {
  s.x -= floor(s.x);
  s.y -= floor(s.y);
  s.z -= floor(s.z);
}
