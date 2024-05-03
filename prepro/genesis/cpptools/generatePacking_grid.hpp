#ifndef GENERATEPACKING_GRID_HPP
#define GENERATEPACKING_GRID_HPP

#include <cmath>
#include <fstream>
#include <iostream>

#include "addParticle.hpp"

#include "message.hpp"
#include "quat.hpp"


// This function generates a grid of particles and writes their data to an output stream. 
// The function takes as input the output stream, a name, the origin and size of the box, 
// the number of particles in each direction, a group number, a homothety factor, and a flag 
// to randomize the orientation of the particles. The function uses nested loops to iterate 
// over each particle in the grid, calculates its position and orientation, and
// writes this data to the output stream. The function finally returns the total number of particles.
//
int generatePacking_grid(std::ostream& os, const char* name, vec3r& origBox, vec3r& boxSize, vec3i& n, int group,
                         double homothety, int randQ = 1) {
  using namespace std;
  int Np = n.x * n.y * n.z;

  quat Qclust;

  int clustID = 0;

  vec3r step(boxSize.x / (double)n.y, boxSize.x / (double)n.y, boxSize.z / (double)n.z);
  vec3r orig = origBox + vec3r(0.5 * step.x, 0.5 * step.y, 0.5 * step.z);

  for (int iy = 0; iy < n.y; iy++) {
    for (int ix = 0; ix < n.x; ix++) {
      for (int iz = 0; iz < n.z; iz++) {

        vec3r gridPos = orig + vec3r(ix * step.x, iy * step.y, iz * step.z);
        if (randQ == 1) {
          Qclust.randomize();
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

#endif /* end of include guard: GENERATEPACKING_GRID_HPP */
