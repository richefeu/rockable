#ifndef GENERATESHAPE_STICK_HPP
#define GENERATESHAPE_STICK_HPP

#include <cmath>
#include <iostream>

void generateShape_stick(std::ostream& os, const char* name, double L, double radius) {
  using namespace std;

  double vol = (4.0 / 3.0) * M_PI * radius * radius * radius + M_PI * radius * radius * L;
  double I_m = radius * radius / 4.0 + (L * L / 12.0);
  double I_m_axis = (1.0 / 2.0) * radius * radius;
  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "preCompDone y" << endl;
  os << "nv 2" << endl;
  os << -L / 2 << " 0 0" << endl;
  os << L / 2 << " 0 0" << endl;
  os << "ne 1" << endl;
  os << "0 1" << endl;
  os << "nf 0" << endl;
  os << "obb.extent " << 0.5 * L + radius << " " << radius << " " << radius << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << vol << endl;
  os << "I/m " << I_m_axis << ' ' << I_m << ' ' << I_m << endl;
  os << ">" << endl;
}


#endif /* end of include guard: GENERATESHAPE_STICK_HPP */
