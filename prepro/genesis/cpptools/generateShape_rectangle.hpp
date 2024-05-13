#ifndef GENERATESHAPE_RECTANGLE_XZ_HPP
#define GENERATESHAPE_RECTANGLE_XZ_HPP

#include <iostream>

#include "vec3.hpp"

void generateShape_rectangle_xz(std::ostream& os, const char* name, double radius, double sidex, double sidez) {
  //using namespace std;
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