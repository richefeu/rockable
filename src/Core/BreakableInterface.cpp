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

#include <algorithm>
#include <cmath>
#include <utility>

#include "BreakableInterface.hpp"
#include "Interaction.hpp"
#include "vec2.hpp"

BreakableInterface::BreakableInterface()
    : i(0), j(0), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}

BreakableInterface::BreakableInterface(size_t I, size_t J)
    : i(I), j(J), kn(0.0), kt(0.0), kr(0.0), fn0(0.0), ft0(0.0), mom0(0.0), power(2.0), dn0(0.0), isInner(1) {}

/**
 * @brief Computes the area of the interface from the current positions of its bonds
 *        (the bonds are supposed to be at the corners of a convex polygon).
 *
 * The bonds are not created in the order of the corners: they are sorted by their angle
 * around their centre, in the plane of the interface (normal to the mean normal of the
 * bonds), and the area is that of the polygon they then make.
 */
void BreakableInterface::computeAreaFromConcernedBonds() {
  area = 0.0;
  size_t n = concernedBonds.size();
  if (n < 3) {
    return;  // Not a valid polygon
  }

  vec3r c, normal;
  for (Interaction* I : concernedBonds) {
    c += I->pos;
    normal += (I->n * concernedBonds[0]->n >= 0.0) ? I->n : -I->n;
  }
  c /= (double)n;
  if (norm2(normal) == 0.0) {
    return;
  }
  normal.normalize();
  vec3r e1 = (std::fabs(normal.x) < 0.9) ? cross(normal, vec3r(1.0, 0.0, 0.0)) : cross(normal, vec3r(0.0, 1.0, 0.0));
  e1.normalize();
  vec3r e2 = cross(normal, e1);

  std::vector<std::pair<double, vec2r>> corners(n);
  for (size_t b = 0; b < n; ++b) {
    vec3r d = concernedBonds[b]->pos - c;
    vec2r p(d * e1, d * e2);
    corners[b] = std::make_pair(atan2(p.y, p.x), p);
  }
  std::sort(corners.begin(), corners.end(),
            [](const std::pair<double, vec2r>& a, const std::pair<double, vec2r>& b) { return a.first < b.first; });

  double twiceArea = 0.0;
  for (size_t b = 0; b < n; ++b) {
    const vec2r& p = corners[b].second;
    const vec2r& q = corners[(b + 1) % n].second;
    twiceArea += p.x * q.y - q.x * p.y;
  }
  area = 0.5 * std::fabs(twiceArea);
}




