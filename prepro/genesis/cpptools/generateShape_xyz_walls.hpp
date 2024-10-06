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

#ifndef GENERATESHAPE_XYZ_WALLS_HPP
#define GENERATESHAPE_XYZ_WALLS_HPP

#include <cmath>
#include <iostream>

#include "vec3.hpp"

void generateShape_xyz_walls(std::ostream& os, double LX, double LY, double LZ, double Rw) {
  using namespace std;

  vec3r ext(0.5 * LX + Rw, 0.5 * LY + Rw, 0.5 * LZ + Rw);
  os << "<" << endl;
  os << "name x-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << "0 " << ext.y << " " << ext.z << endl;
  os << "0 " << ext.y << " " << -ext.z << endl;
  os << "0 " << -ext.y << " " << -ext.z << endl;
  os << "0 " << -ext.y << " " << ext.z << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << Rw << " " << ext.y + Rw << " " << ext.z + Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LX + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;

  os << "<" << endl;
  os << "name y-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << ext.x << " 0 " << ext.z << endl;
  os << ext.x << " 0 " << -ext.z << endl;
  os << -ext.x << " 0 " << -ext.z << endl;
  os << -ext.x << " 0 " << ext.z << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << ext.x + Rw << " " << Rw << " " << ext.z + Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LY + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;

  os << "<" << endl;
  os << "name z-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << ext.x << " " << ext.y << " 0" << endl;
  os << ext.x << " " << -ext.y << " 0" << endl;
  os << -ext.x << " " << -ext.y << " 0" << endl;
  os << -ext.x << " " << ext.y << " 0" << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << ext.x + Rw << " " << ext.y + Rw << " " << Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LZ + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;
}

void generateShape_xyz_walls(std::ostream& os, vec3r & size, double Rw) {
  generateShape_xyz_walls(os, size.x, size.y, size.z, Rw);
}


#endif /* end of include guard: GENERATESHAPE_XYZ_WALLS_HPP */
