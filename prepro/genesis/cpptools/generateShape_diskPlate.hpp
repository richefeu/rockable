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
