#ifndef GENERATESHAPE_CUBOID_HPP
#define GENERATESHAPE_CUBOID_HPP

#include <iostream>

#include "vec3.hpp"

// sideSize = external sizes of the cuboid block (it includes the radius)
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

  os << "<" << '\n';
  os << "name " << name << '\n';
  os << "radius " << radius << '\n';
  os << "preCompDone y" << '\n';
  os << "nv 8" << '\n';
  os << lx << " " << ly << " " << -lz << '\n';
  os << -lx << " " << ly << " " << -lz << '\n';
  os << -lx << " " << -ly << " " << -lz << '\n';
  os << lx << " " << -ly << " " << -lz << '\n';
  os << lx << " " << ly << " " << lz << '\n';
  os << -lx << " " << ly << " " << lz << '\n';
  os << -lx << " " << -ly << " " << lz << '\n';
  os << lx << " " << -ly << " " << lz << '\n';
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
  os << "obb.extent " << 0.5 * sideSize.x << ' ' << 0.5 * sideSize.y << ' ' << 0.5 * sideSize.z << '\n';
  os << "obb.e1 1 0 0" << '\n';
  os << "obb.e2 0 1 0" << '\n';
  os << "obb.e3 0 0 1" << '\n';
  os << "obb.center 0 0 0" << '\n';
  os << "volume " << vol << '\n';
  os << "I/m " << Ix_m << ' ' << Iy_m << ' ' << Iz_m << '\n';
  os << ">" << '\n';
  os << std::flush;
}

#endif /* end of include guard: GENERATESHAPE_CUBOID_HPP */