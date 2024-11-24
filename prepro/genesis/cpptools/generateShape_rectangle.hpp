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

#ifndef GENERATESHAPE_RECTANGLE_XZ_HPP
#define GENERATESHAPE_RECTANGLE_XZ_HPP

#include <iostream>

#include "vec3.hpp"

void generateShape_rectangle_xy(std::ostream& os, const char* name, double radius, double sidex, double sidey) {
  double lx = 0.5 * sidex - radius;
  double ly = 0.5 * sidey - radius;

  os << "<\n";
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y\n";
  os << "nv 4" << '\n';
  os << lx << " " << ly << " " << 0 << '\n';
  os << lx << " " << -ly << " " << 0 << '\n';
  os << -lx << " " << -ly << " " << 0 << '\n';
  os << -lx << " " << ly << " " << 0 << '\n';
  os << "ne 4" << '\n';
  os << "0 1" << '\n';
  os << "1 2" << '\n';
  os << "2 3" << '\n';
  os << "3 0" << '\n';
  os << "nf 1" << '\n';
  os << "4 0 1 2 3" << '\n';
  os << "obb.extent " << 0.5 * sidex << ' ' << 0.5 * sidey  << ' ' << radius << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << sidex * sidey * 2.0 * radius << '\n'; // bidon
  os << "I/m " << 1.0 << ' ' << 1.0 << ' ' << 1.0 << '\n'; // bidon
  os << ">\n";
  os << std::flush;
}

void generateShape_rectangle_xz(std::ostream& os, const char* name, double radius, double sidex, double sidez) {
  double lx = 0.5 * sidex - radius;
  double lz = 0.5 * sidez - radius;

  os << "<\n";
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y\n";
  os << "nv 4" << '\n';
  os << lx << " " << 0 << " " << lz << '\n';
  os << lx << " " << 0 << " " << -lz << '\n';
  os << -lx << " " << 0 << " " << -lz << '\n';
  os << -lx << " " << 0 << " " << lz << '\n';
  os << "ne 4" << '\n';
  os << "0 1" << '\n';
  os << "1 2" << '\n';
  os << "2 3" << '\n';
  os << "3 0" << '\n';
  os << "nf 1" << '\n';
  os << "4 0 1 2 3" << '\n';
  os << "obb.extent " << 0.5 * sidex << ' ' << radius << ' ' << 0.5 * sidez << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << sidex * sidez * 2.0 * radius << '\n'; // bidon
  os << "I/m " << 1. << ' ' << 1. << ' ' << 1. << '\n'; // bidon
  os << ">\n";
  os << std::flush;
}

#endif /* end of include guard: GENERATESHAPE_RECTANGLE_XZ_HPP */