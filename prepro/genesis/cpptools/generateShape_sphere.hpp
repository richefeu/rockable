#ifndef GENERATESHAPE_SPHERE_HPP
#define GENERATESHAPE_SPHERE_HPP

#include <iostream>
#include <cmath>

void generateShape_sphere(std::ostream& os, const char* name, double radius) {
  using namespace std;

  double vol = (4.0 / 3.0) * M_PI * radius * radius * radius;
  double I_m = (2.0 / 5.0) * radius * radius;
  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "preCompDone y" << endl;
  os << "nv 1" << endl;
  os << "0 0 0" << endl;
  os << "ne 0" << endl;
  os << "nf 0" << endl;
  os << "obb.extent 1 1 1" << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << vol << endl;
  os << "I/m " << I_m << ' ' << I_m << ' ' << I_m << endl;
  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_SPHERE_HPP */
