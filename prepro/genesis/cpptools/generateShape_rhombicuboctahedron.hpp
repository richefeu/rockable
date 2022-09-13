#ifndef GENERATESHAPE_RHOMBICUBOCTAHEDRON_HPP
#define GENERATESHAPE_RHOMBICUBOCTAHEDRON_HPP

#include <iostream>

#include "vec3.hpp"

// sideSize = external sizes of the rhombicuboctahedron block (it includes the radius)
void generateShape_rhombicuboctahedron(std::ostream& os, const char* name, double radius, vec3r& sideSize) {
  using namespace std;
  double lx = 0.5 * sideSize.x - radius;
  double ly = 0.5 * sideSize.y - radius;
  double lz = 0.5 * sideSize.z - radius;
  double ix = lx / (1.0 + sqrt(2.0));
  double iy = ly / (1.0 + sqrt(2.0));
  double iz = lz / (1.0 + sqrt(2.0));

  /*
  double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  voidVol += (M_PI * radius * radius * (sideSize.x - 2.0 * radius));
  voidVol += (M_PI * radius * radius * (sideSize.y - 2.0 * radius));
  voidVol += (M_PI * radius * radius * (sideSize.z - 2.0 * radius));
  double vol = sideSize.x * sideSize.y * sideSize.z - voidVol;

  double Ix_m = (sideSize.y * sideSize.y + sideSize.z * sideSize.z) / 12.0;
  double Iy_m = (sideSize.z * sideSize.z + sideSize.x * sideSize.x) / 12.0;
  double Iz_m = (sideSize.x * sideSize.x + sideSize.y * sideSize.y) / 12.0;
  */

  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "preCompDone n" << endl;
  os << "MCnstep 100000" << endl;

  os << "nv 24" << endl;
  os << lx << " " << -iy << " " << iz << endl;
  os << lx << " " << -iy << " " << -iz << endl;
  os << lx << " " << iy << " " << -iz << endl;
  os << lx << " " << iy << " " << iz << endl;

  os << ix << " " << -ly << " " << iz << endl;
  os << ix << " " << -ly << " " << -iz << endl;
  os << ix << " " << -iy << " " << -lz << endl;
  os << ix << " " << iy << " " << -lz << endl;
  os << ix << " " << ly << " " << -iz << endl;
  os << ix << " " << ly << " " << iz << endl;
  os << ix << " " << iy << " " << lz << endl;
  os << ix << " " << -iy << " " << lz << endl;

  os << -ix << " " << -ly << " " << iz << endl;
  os << -ix << " " << -ly << " " << -iz << endl;
  os << -ix << " " << -iy << " " << -lz << endl;
  os << -ix << " " << iy << " " << -lz << endl;
  os << -ix << " " << ly << " " << -iz << endl;
  os << -ix << " " << ly << " " << iz << endl;
  os << -ix << " " << iy << " " << lz << endl;
  os << -ix << " " << -iy << " " << lz << endl;

  os << -lx << " " << -iy << " " << iz << endl;
  os << -lx << " " << -iy << " " << -iz << endl;
  os << -lx << " " << iy << " " << -iz << endl;
  os << -lx << " " << iy << " " << iz << endl;

  os << "ne 48" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;

  os << "4 5" << endl;
  os << "5 6" << endl;
  os << "6 7" << endl;
  os << "7 8" << endl;
  os << "8 9" << endl;
  os << "9 10" << endl;
  os << "10 11" << endl;
  os << "11 4" << endl;

  os << "12 13" << endl;
  os << "13 14" << endl;
  os << "14 15" << endl;
  os << "15 16" << endl;
  os << "16 17" << endl;
  os << "17 18" << endl;
  os << "18 19" << endl;
  os << "19 12" << endl;

  os << "20 21" << endl;
  os << "21 22" << endl;
  os << "22 23" << endl;
  os << "23 20" << endl;

  os << "0 11" << endl;
  os << "0 4" << endl;
  os << "1 5" << endl;
  os << "1 6" << endl;
  os << "2 7" << endl;
  os << "2 8" << endl;
  os << "3 9" << endl;
  os << "3 10" << endl;

  os << "4 12" << endl;
  os << "5 13" << endl;
  os << "6 14" << endl;
  os << "7 15" << endl;
  os << "8 16" << endl;
  os << "9 17" << endl;
  os << "10 18" << endl;
  os << "11 19" << endl;

  os << "20 19" << endl;
  os << "20 12" << endl;
  os << "21 13" << endl;
  os << "21 14" << endl;
  os << "22 15" << endl;
  os << "22 16" << endl;
  os << "23 17" << endl;
  os << "23 18" << endl;

  os << "nf 26" << endl;

  os << "4 0 1 2 3" << endl;
  os << "4 4 5 13 12" << endl;
  os << "4 6 7 15 14" << endl;
  os << "4 8 16 17 9" << endl;
  os << "4 11 10 18 19" << endl;
  os << "4 20 21 22 23" << endl;

  os << "4 0 1 5 4" << endl;
  os << "4 1 2 7 6" << endl;
  os << "4 2 8 9 3" << endl;
  os << "4 0 3 10 11" << endl;

  os << "4 5 6 14 13" << endl;
  os << "4 7 8 16 15" << endl;
  os << "4 9 10 18 17" << endl;
  os << "4 4 11 19 12" << endl;

  os << "4 12 13 21 20" << endl;
  os << "4 14 15 22 21" << endl;
  os << "4 16 17 23 22" << endl;
  os << "4 18 19 20 23" << endl;

  os << "3 0 11 4" << endl;
  os << "3 1 6 5" << endl;
  os << "3 2 7 8" << endl;
  os << "3 3 9 10" << endl;

  os << "3 20 12 19" << endl;
  os << "3 21 13 14" << endl;
  os << "3 22 15 16" << endl;
  os << "3 23 17 18" << endl;

  /*
  os << "obb.extent " << 0.5 * sideSize.x << ' ' << 0.5 * sideSize.y << ' ' << 0.5 * sideSize.z << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << vol << endl;
  os << "I/m " << Ix_m << ' ' << Iy_m << ' ' << Iz_m << endl;
  */
  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_RHOMBICUBOCTAHEDRON_HPP */