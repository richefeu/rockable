#ifndef GENERATESHAPE_XYZ_WALLS_HPP
#define GENERATESHAPE_XYZ_WALLS_HPP

#include <cmath>
#include <iostream>

#include "vec3.hpp"

void generateShape_xyz_walls(std::ostream& os, double LX, double LY, double LZ, double Rw) {
  using namespace std;

  vec3r ext(0.5 * LX + Rw, 0.5 * LY + Rw, 0.5 * LZ + Rw);
  os << "<" << endl;
  os << "name x-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << "0 " << ext.y << " " << ext.z << endl;
  os << "0 " << ext.y << " " << -ext.z << endl;
  os << "0 " << -ext.y << " " << -ext.z << endl;
  os << "0 " << -ext.y << " " << ext.z << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << Rw << " " << ext.y + Rw << " " << ext.z + Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LX + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;

  os << "<" << endl;
  os << "name y-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << ext.x << " 0 " << ext.z << endl;
  os << ext.x << " 0 " << -ext.z << endl;
  os << -ext.x << " 0 " << -ext.z << endl;
  os << -ext.x << " 0 " << ext.z << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << ext.x + Rw << " " << Rw << " " << ext.z + Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LY + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;

  os << "<" << endl;
  os << "name z-wall" << endl;
  os << "radius " << Rw << endl;
  os << "preCompDone y" << endl;
  os << "nv 4" << endl;
  os << ext.x << " " << ext.y << " 0" << endl;
  os << ext.x << " " << -ext.y << " 0" << endl;
  os << -ext.x << " " << -ext.y << " 0" << endl;
  os << -ext.x << " " << ext.y << " 0" << endl;
  os << "ne 4" << endl;
  os << "0 1" << endl;
  os << "1 2" << endl;
  os << "2 3" << endl;
  os << "3 0" << endl;
  os << "nf 1" << endl;
  os << "4 0 1 2 3" << endl;
  os << "obb.extent " << ext.x + Rw << " " << ext.y + Rw << " " << Rw << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 0 0" << endl;
  os << "volume " << (LZ + 2 * Rw) * 2 * Rw << endl;
  os << "I/m 1.0 1.0 1.0" << endl;
  os << ">" << endl;
}

void generateShape_xyz_walls(std::ostream& os, vec3r & size, double Rw) {
  generateShape_xyz_walls(os, size.x, size.y, size.z, Rw);
}


#endif /* end of include guard: GENERATESHAPE_XYZ_WALLS_HPP */
