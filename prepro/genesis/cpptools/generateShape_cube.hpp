#ifndef GENERATESHAPE_CUBE_HPP
#define GENERATESHAPE_CUBE_HPP

#include <iostream>

// sideSize = external size of the cube (it includes the radius)
void generateShape_cube(std::ostream& os, const char* name, double radius, double sideSize) {
  using namespace std;
  double l = 0.5 * sideSize - radius;
  double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  voidVol += 3.0 * (M_PI * radius * radius * (sideSize - 2.0 * radius));
  double vol = sideSize * sideSize * sideSize - voidVol;
  double I_m = (sideSize * sideSize) / 6.0;

  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "preCompDone y" << endl;
  os << "nv 8" << endl;
  os << l << " " << l << " " << -l << endl;
  os << -l << " " << l << " " << -l << endl;
  os << -l << " " << -l << " " << -l << endl;
  os << l << " " << -l << " " << -l << endl;
  os << l << " " << l << " " << l << endl;
  os << -l << " " << l << " " << l << endl;
  os << -l << " " << -l << " " << l << endl;
  os << l << " " << -l << " " << l << endl;
  os << "ne 12" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "4 5" << endl;
  os << "5 6" << endl;
  os << "6 7" << endl;
  os << "7 4" << endl;
  os << "0 4" << endl;
  os << "1 5" << endl;
  os << "2 6" << endl;
  os << "3 7" << endl;
  os << "nf 6" << endl;
  os << "4 0 1 2 3" << endl;
  os << "4 4 5 6 7" << endl;
  os << "4 0 1 5 4" << endl;
  os << "4 2 3 7 6" << endl;
  os << "4 1 2 6 5" << endl;
  os << "4 0 4 7 3" << endl;
  os << "obb.extent " << 0.5 * sideSize << ' ' << 0.5 * sideSize << ' ' << 0.5 * sideSize << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << vol << endl;
  os << "I/m " << I_m << ' ' << I_m << ' ' << I_m << endl;
  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_CUBE_HPP */