#ifndef GENERATESHAPE_CYLINDRICALMOLD_HPP
#define GENERATESHAPE_CYLINDRICALMOLD_HPP

#include <cmath>
#include <iostream>

// radius = Minskowski radius
// height = internal height
// Rin = internal radius
// The origin coordinate is at the bottom inside
// I/m is fake
void generateShape_cylindricalMold(std::ostream& os, const char* name, double radius, double height, int nbSectors,
                                   double Rin) {
  using namespace std;

  double angle = 2.0 * M_PI / nbSectors;

  os << "<" << endl;
  os << "name " << name << endl;
  os << "radius " << radius << endl;
  os << "isSurface" << endl;
  os << "preCompDone y" << endl;
  os << "nv " << 2 * nbSectors << endl;

  // VERTICES
  double R = Rin + radius;
  for (int i = 0; i < nbSectors; i++) {
    double z = R * cos(i * angle);
    double x = R * sin(i * angle);
    os << x << " " << height + radius << " " << z << endl;
  }
  for (int i = 0; i < nbSectors; i++) {
    double z = R * cos(i * angle);
    double x = R * sin(i * angle);
    os << x << " " << -radius << " " << z << endl;
  }
  os << endl;

  // EDGES
  os << "ne " << 3 * nbSectors << endl;
  for (int i = 0; i < nbSectors; i++) {
    int j = i + 1;
    if (j == nbSectors) j = 0;
    os << i << " " << j << endl;
  }
  for (int i = 0; i < nbSectors; i++) {
    int j = nbSectors + i + 1;
    if (j == 2 * nbSectors) j = nbSectors;
    os << nbSectors + i << " " << j << endl;
  }

  for (int i = 0; i < nbSectors; i++) {
    int j = nbSectors + i;
    os << i << " " << j << endl;
  }
  os << endl;

  // FACES
  os << "nf " << nbSectors + 1 << endl;

  for (int n = 0; n < nbSectors; n++) {
    int i = n;
    int j = i + 1;
    if (j == nbSectors) j = 0;
    int k = i + nbSectors + 1;
    if (k == 2 * nbSectors) k = nbSectors;
    int l = i + nbSectors;

    os << "4 " << i << " " << j << " " << k << " " << l << endl;
  }

  os << nbSectors;
  for (int n = 0; n < nbSectors; n++) {
    os << " " << nbSectors + n;
  }
  os << endl;

  os << "obb.extent " << Rin + 2 * radius << " " << 0.5 * height + 2 * radius << " " << Rin + 2 * radius << endl;
  os << "obb.e1 1 0 0" << endl;
  os << "obb.e2 0 1 0" << endl;
  os << "obb.e3 0 0 1" << endl;
  os << "obb.center 0 " << height * 0.5 << " 0" << endl;
  os << "volume " << (2 * Rin * radius + radius * radius) * M_PI * height << endl;
  os << "I/m 1.0 1.0 1.0" << endl;

  os << ">" << endl;
}

#endif /* end of include guard: GENERATESHAPE_CYLINDRICALMOLD_HPP */
