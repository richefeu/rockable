//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

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
