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

#ifndef GENERATEPACKING_GRID_HPP
#define GENERATEPACKING_GRID_HPP

#include <cmath>
#include <fstream>
#include <iostream>

#include "addParticle.hpp"

#include "message.hpp"
#include "quat.hpp"
#include "transformation.hpp"

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp 

// This function generates a grid of particles and writes their data to an output stream.
// The cluter number is incremented, starting with clustID
//
int generatePacking_grid(std::ostream& os, const char* name, vec3r& origBox, vec3r& boxSize, vec3i& n, int group,
                         int clustID, double homothety, int randQ = 1) {
  using namespace std;
  int Np = n.x * n.y * n.z;

  quat Qclust;

  vec3r step(boxSize.x / (double)n.x, boxSize.y / (double)n.y, boxSize.z / (double)n.z);
  vec3r orig = origBox + vec3r(0.5 * step.x, 0.5 * step.y, 0.5 * step.z);

  for (int iy = 0; iy < n.y; iy++) {
    for (int ix = 0; ix < n.x; ix++) {
      for (int iz = 0; iz < n.z; iz++) {

        vec3r gridPos = orig + vec3r(ix * step.x, iy * step.y, iz * step.z);
        globalTransformation.apply(gridPos);
        if (randQ == 1) {
          Qclust.randomize();
        } else {
          Qclust = individualParticleRotation;
        }

        msg::bestPrecision(os);
        os << name << ' ' << group << ' ' << clustID << ' ' << homothety << ' ' << gridPos << "  0 0 0  0 0 0  "
           << Qclust << "  0 0 0  0 0 0\n";

        clustID++;
        msg::normalPrecision(os);
      }
    }
  }

  return Np;
}

// This function generates a grid of particles and writes their data to an output stream.
// The cluter number is clustID for all generated particle
//
int generatePacking_grid_clust(std::ostream& os, const char* name, vec3r& origBox, vec3r& boxSize, vec3i& n, int group,
                         int clustID, double homothety, int randQ = 1) {
  using namespace std;
  int Np = n.x * n.y * n.z;

  quat Qclust;

  vec3r step(boxSize.x / (double)n.x, boxSize.y / (double)n.y, boxSize.z / (double)n.z);
  vec3r orig = origBox + vec3r(0.5 * step.x, 0.5 * step.y, 0.5 * step.z);

  for (int iy = 0; iy < n.y; iy++) {
    for (int ix = 0; ix < n.x; ix++) {
      for (int iz = 0; iz < n.z; iz++) {

        vec3r gridPos = orig + vec3r(ix * step.x, iy * step.y, iz * step.z);
        globalTransformation.apply(gridPos);
        if (randQ == 1) {
          Qclust.randomize();
        } else {
          Qclust = individualParticleRotation;
        }

        msg::bestPrecision(os);
        os << name << ' ' << group << ' ' << clustID << ' ' << homothety << ' ' << gridPos << "  0 0 0  0 0 0  "
           << Qclust << "  0 0 0  0 0 0\n";

        msg::normalPrecision(os);
      }
    }
  }

  return Np;
}

#endif /* end of include guard: GENERATEPACKING_GRID_HPP */
