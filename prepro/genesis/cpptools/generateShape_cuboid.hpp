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

#ifndef GENERATESHAPE_CUBOID_HPP
#define GENERATESHAPE_CUBOID_HPP

#include <iostream>

#include "vec3.hpp"

// sideSize = external sizes of the cuboid block (it includes the radius)
void generateShape_cuboid(std::ostream& os, const char* name, double radius, vec3r& sideSize) {
  //using namespace std;
  double lx = 0.5 * sideSize.x - radius;
	double ly = 0.5 * sideSize.y - radius;
	double lz = 0.5 * sideSize.z - radius;
  double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  voidVol += (M_PI * radius * radius * (sideSize.x - 2.0 * radius));
	voidVol += (M_PI * radius * radius * (sideSize.y - 2.0 * radius));
	voidVol += (M_PI * radius * radius * (sideSize.z - 2.0 * radius));
  double vol = sideSize.x * sideSize.y * sideSize.z - voidVol;
  double Ix_m = (sideSize.y * sideSize.y + sideSize.z * sideSize.z) / 12.0;
	double Iy_m = (sideSize.z * sideSize.z + sideSize.x * sideSize.x) / 12.0;
	double Iz_m = (sideSize.x * sideSize.x + sideSize.y * sideSize.y) / 12.0;

  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y" << '\n';
  os << "nv 8" << '\n';
  os << lx << " " << ly << " " << -lz << '\n';
  os << -lx << " " << ly << " " << -lz << '\n';
  os << -lx << " " << -ly << " " << -lz << '\n';
  os << lx << " " << -ly << " " << -lz << '\n';
  os << lx << " " << ly << " " << lz << '\n';
  os << -lx << " " << ly << " " << lz << '\n';
  os << -lx << " " << -ly << " " << lz << '\n';
  os << lx << " " << -ly << " " << lz << '\n';
  os << "ne 12" << '\n';
  os << "0 1" << '\n';
  os << "1 2" << '\n';
  os << "2 3" << '\n';
  os << "3 0" << '\n';
  os << "4 5" << '\n';
  os << "5 6" << '\n';
  os << "6 7" << '\n';
  os << "7 4" << '\n';
  os << "0 4" << '\n';
  os << "1 5" << '\n';
  os << "2 6" << '\n';
  os << "3 7" << '\n';
  os << "nf 6" << '\n';
  os << "4 0 1 2 3" << '\n';
  os << "4 4 5 6 7" << '\n';
  os << "4 0 1 5 4" << '\n';
  os << "4 2 3 7 6" << '\n';
  os << "4 1 2 6 5" << '\n';
  os << "4 0 4 7 3" << '\n';
  os << "obb.extent " << 0.5 * sideSize.x << ' ' << 0.5 * sideSize.y << ' ' << 0.5 * sideSize.z << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << vol << '\n';
  os << "I/m " << Ix_m << ' ' << Iy_m << ' ' << Iz_m << '\n';
  os << ">" << '\n';
  os << std::flush;
}

#endif /* end of include guard: GENERATESHAPE_CUBOID_HPP */