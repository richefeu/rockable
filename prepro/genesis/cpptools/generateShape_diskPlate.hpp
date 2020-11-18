#ifndef GENERATESHAPE_DISKPLATE_HPP
#define GENERATESHAPE_DISKPLATE_HPP

#include <cmath>
#include <iostream>

// radius = Minskowski radius
// Rout = external radius
// The origin coordinate is in the middle
// I/m is a fake one
void generateShape_diskPlate(std::ostream& os, const char* name, double radius, int nbSectors, double Rout) {
  using namespace std;

  double angle = 2.0 * M_PI / (double)nbSectors;

  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "isSurface" << endl;
  os << "preCompDone y" << endl;
  os << "nv " << nbSectors << endl;

  // VERTICES
  double R = Rout - radius;
  for (int i = 0; i < nbSectors; i++) {
    double z = R * cos((double)i * angle);
    double x = R * sin((double)i * angle);
    os << x << " " << 0.0 << " " << z << endl;
  }
  os << endl;

  // EDGES
  os << "ne " << nbSectors << endl;
  for (int i = 0; i < nbSectors; i++) {
    int j = i + 1;
    if (j == nbSectors) j = 0;
    os << i << " " << j << endl;
  }
  os << endl;

  // FACE
  os << "nf " << 1 << endl;
  os << nbSectors;
  for (int n = 0; n < nbSectors; n++) {
    os << " " << n;
  }
  os << endl;

  os << "obb.extent " << Rout << " " << radius << " " << Rout << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 " << 0 << " 0" << endl;
  os << "volume " << M_PI * Rout * Rout * (2 * radius) << endl;
  double I_m = Rout * Rout / 4.0 + (4 * radius * radius / 12.0);
  double I_m_axis = (1.0 / 2.0) * Rout * Rout;
  os << "I/m " << I_m << " " << I_m_axis << " " << I_m << endl;

  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_DISKPLATE_HPP */
