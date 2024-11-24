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

#ifndef GENERATEPACKING_WALLBOX_HPP
#define GENERATEPACKING_WALLBOX_HPP

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp 

// TODO: for now the transformations are not applied

#include <iostream>

// Lengths LX, LY and LZ are inside (walls are placed outside this cube)
// The rank-order is xmin, xmax, ymin, ymax, zmin, zmax
// TODO: add an origin
int generatePacking_wallBox(std::ostream& os, int group, double LX, double LY, double LZ, double Rw) {

  // normal X
  os << "x-wall " << group << " 0 1 " << -Rw << " " << 0.5 * LY << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;
  os << "x-wall " << group << " 0 1 " << LX + Rw << " " << 0.5 * LY << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;

  // normal Y
  os << "y-wall " << group << " 0 1 " << 0.5 * LX << " " << -Rw << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;
  os << "y-wall " << group << " 0 1 " << 0.5 * LX << " " << LY + Rw << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;

  // normal Z
  os << "z-wall " << group << " 0 1 " << 0.5 * LX << " " << 0.5 * LY << " " << -Rw
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;
  os << "z-wall " << group << " 0 1 " << 0.5 * LX << " " << 0.5 * LY << " " << LZ + Rw
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << std::endl;

  return 6;
}

#endif /* end of include guard: GENERATEPACKING_WALLBOX_HPP */
