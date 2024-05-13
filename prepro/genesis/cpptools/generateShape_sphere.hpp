#ifndef GENERATESHAPE_SPHERE_HPP
#define GENERATESHAPE_SPHERE_HPP

#include <iostream>
#include <cmath>

void generateShape_sphere(std::ostream& os, const char* name, double radius) {
  double vol = (4.0 / 3.0) * M_PI * radius * radius * radius;
  double I_m = (2.0 / 5.0) * radius * radius;
  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y" << '\n';
  os << "nv 1" << '\n';
  os << "0 0 0" << '\n';
  os << "ne 0" << '\n';
  os << "nf 0" << '\n';
  os << "obb.extent 1 1 1" << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << vol << '\n';
  os << "I/m " << I_m << ' ' << I_m << ' ' << I_m << '\n';
  os << ">" << std::endl;
}

#endif /* end of include guard: GENERATESHAPE_SPHERE_HPP */
