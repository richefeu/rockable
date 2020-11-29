#ifndef GENERATESHAPE_CUBOID_HPP
#define GENERATESHAPE_CUBOID_HPP

#include <iostream>

#include "vec3.hpp"

// sideSize = external size of the cuboid block (it includes the radius)
void generateShape_cuboid(std::ostream& os, const char* name, double radius, vec3r& sideSize) {
  using namespace std;
  double lx = 0.5 * sideSize.x - radius;
	double ly = 0.5 * sideSize.y - radius;
	double lz = 0.5 * sideSize.z - radius;
  double voidVol = (8.0 - (4.0 / 3.0) * M_PI) * radius * radius * radius;
  voidVol += (M_PI * radius * radius * (sideSize.x - 2.0 * radius));
	voidVol += (M_PI * radius * radius * (sideSize.y - 2.0 * radius));
	voidVol += (M_PI * radius * radius * (sideSize.z - 2.0 * radius));
  double vol = sideSize.x * sideSize.y * sideSize.z - voidVol;
  double Ix_m = (sideSize.y * sideSize.y + sideSize.z * sideSize.z) / 12.0;
	double Iy_m = (sideSize.z * sideSize.z + sideSize.x * sideSize.x) / 12.0;
	double Iz_m = (sideSize.x * sideSize.x + sideSize.y * sideSize.y) / 12.0;

  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "preCompDone y" << endl;
  os << "nv 8" << endl;
  os << lx << " " << ly << " " << -lz << endl;
  os << -lx << " " << ly << " " << -lz << endl;
  os << -lx << " " << -ly << " " << -lz << endl;
  os << lx << " " << -ly << " " << -lz << endl;
  os << lx << " " << ly << " " << lz << endl;
  os << -lx << " " << ly << " " << lz << endl;
  os << -lx << " " << -ly << " " << lz << endl;
  os << lx << " " << -ly << " " << lz << endl;
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
  os << "obb.extent " << 0.5 * sideSize.x << ' ' << 0.5 * sideSize.y << ' ' << 0.5 * sideSize.z << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << vol << endl;
  os << "I/m " << Ix_m << ' ' << Iy_m << ' ' << Iz_m << endl;
  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_CUBOID_HPP */