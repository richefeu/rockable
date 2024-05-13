#ifndef GENERATESHAPE_CUBE_HPP
#define GENERATESHAPE_CUBE_HPP

#include <iostream>

// sideSize = external size of the cube (it includes the radius)
void generateShape_cube(std::ostream& os, const char* name, double radius, double sideSize) {

  double l = 0.5 * sideSize - radius;
  double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  voidVol += 3.0 * (M_PI * radius * radius * (sideSize - 2.0 * radius));
  double vol = sideSize * sideSize * sideSize - voidVol;
  double I_m = (sideSize * sideSize) / 6.0;

  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y" << '\n';
  os << "nv 8" << '\n';
  os << l << " " << l << " " << -l << '\n';
  os << -l << " " << l << " " << -l << '\n';
  os << -l << " " << -l << " " << -l << '\n';
  os << l << " " << -l << " " << -l << '\n';
  os << l << " " << l << " " << l << '\n';
  os << -l << " " << l << " " << l << '\n';
  os << -l << " " << -l << " " << l << '\n';
  os << l << " " << -l << " " << l << '\n';
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
  os << "obb.extent " << 0.5 * sideSize << ' ' << 0.5 * sideSize << ' ' << 0.5 * sideSize << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << vol << '\n';
  os << "I/m " << I_m << ' ' << I_m << ' ' << I_m << '\n';
  os << ">" << std::endl;
}

#endif /* end of include guard: GENERATESHAPE_CUBE_HPP */