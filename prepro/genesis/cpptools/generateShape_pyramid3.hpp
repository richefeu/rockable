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

#ifndef GENERATESHAPE_PYRAMID3_HPP
#define GENERATESHAPE_PYRAMID3_HPP

#include <cmath>
#include <iostream>

// sideSize = external size of the cube (it includes the radius)
void generateShape_pyramid3(std::ostream& os, const char* name, double radius, double sideSize) {

  double l = 0.5 * sideSize - radius;
  double sqrt2 = sqrt(2.0);
  // double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  // voidVol += 3.0 * (M_PI * radius * radius * (sideSize - 2.0 * radius));
  // double vol = 0.118 * sideSize * sideSize * sideSize ;
  // double I_m = (sideSize * sideSize) / 6.0;

  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone n" << '\n';
  os << "MCnstep 100000" << '\n';
  os << "nv 4" << '\n';
  os << -l << " " << 0 << " " << -0.5 * sqrt2 * l << '\n';
  os << l << " " << 0 << " " << -0.5 * sqrt2 * l << '\n';
  os << 0 << " " << -l << " " << 0.5 * sqrt2 * l << '\n';
  os << 0 << " " << l << " " << 0.5 * sqrt2 * l << '\n';
  os << "ne 6" << '\n';
  os << "0 1" << '\n';
  os << "0 2" << '\n';
  os << "0 3" << '\n';
  os << "1 2" << '\n';
  os << "1 3" << '\n';
  os << "2 3" << '\n';
  os << "nf 4" << '\n';
  os << "3 0 1 2" << '\n';
  os << "3 0 2 3" << '\n';
  os << "3 2 1 3" << '\n';
  os << "3 1 0 3" << '\n';
  os << "obb.extent " << l + radius << ' ' << l + radius << ' ' << 0.5 * sqrt2 * l + radius << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << 1.0 << '\n';
  os << "I/m " << 1.0 << ' ' << 1.0 << ' ' << 1.0 << '\n';
  os << ">" << std::endl;
}

#endif /* end of include guard: GENERATESHAPE_PYRAMID3_HPP */
