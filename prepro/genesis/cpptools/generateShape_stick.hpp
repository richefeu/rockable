#ifndef GENERATESHAPE_STICK_HPP
#define GENERATESHAPE_STICK_HPP

#include <cmath>
#include <iostream>

void generateShape_stick(std::ostream& os, const char* name, double L, double radius) {

  double vol = (4.0 / 3.0) * M_PI * radius * radius * radius + M_PI * radius * radius * L;
  double I_m = radius * radius / 4.0 + (L * L / 12.0);
  double I_m_axis = (1.0 / 2.0) * radius * radius;
  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y" << '\n';
  os << "nv 2" << '\n';
  os << -L / 2 << " 0 0" << '\n';
  os << L / 2 << " 0 0" << '\n';
  os << "ne 1" << '\n';
  os << "0 1" << '\n';
  os << "nf 0" << '\n';
  os << "obb.extent " << 0.5 * L + radius << " " << radius << " " << radius << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << vol << '\n';
  os << "I/m " << I_m_axis << ' ' << I_m << ' ' << I_m << '\n';
  os << ">" << std::endl;
}


#endif /* end of include guard: GENERATESHAPE_STICK_HPP */
