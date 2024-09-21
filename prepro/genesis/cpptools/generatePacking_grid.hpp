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
