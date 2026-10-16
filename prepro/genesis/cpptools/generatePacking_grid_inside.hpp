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

#ifndef GENERATEPACKING_GRID_INSIDE_HPP
#define GENERATEPACKING_GRID_INSIDE_HPP

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "message.hpp"
#include "quat.hpp"
#include "transformation.hpp"
#include "vec3.hpp"

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp

// A grid of particles cut by a shape (a 'cookie cutter'): only the positions of
// the grid that are inside the shape are kept. The shape is read in a shape
// file; it is placed at the origin, in the frame of the file, scaled by a
// homothety. Its faces must form a closed surface.

struct GridCutter {
  std::string name;
  std::vector<vec3r> vertex;
  std::vector<std::vector<size_t>> face;
};

// Reads the vertices and the faces of the shape 'name' in the shape file 'fileName'
bool readGridCutter(const char* fileName, const std::string& name, GridCutter& C) {
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::cerr << msg::warn() << "@readGridCutter, cannot read " << fileName << msg::normal() << '\n';
    return false;
  }
  std::string line;
  bool inShape = false;
  bool found = false;
  while (std::getline(file, line)) {
    std::istringstream is(line);
    std::string key;
    if (!(is >> key)) continue;
    if (key == "<") {
      inShape = false;
    } else if (key == ">") {
      if (inShape) {
        found = true;
        break;
      }
    } else if (key == "name") {
      std::string n;
      is >> n;
      inShape = (n == name);
      if (inShape) {
        C.vertex.clear();
        C.face.clear();
      }
    } else if (key == "nv" || key == "ne" || key == "nf") {
      size_t n = 0;
      is >> n;
      for (size_t i = 0; i < n; i++) {
        std::getline(file, line);
        if (!inShape) continue;
        std::istringstream il(line);
        if (key == "nv") {
          vec3r v;
          il >> v;
          C.vertex.push_back(v);
        } else if (key == "nf") {
          size_t nb = 0;
          il >> nb;
          std::vector<size_t> f(nb);
          for (size_t k = 0; k < nb; k++) il >> f[k];
          C.face.push_back(f);
        }
      }
    }
  }
  if (!found || C.face.empty()) {
    std::cerr << msg::warn() << "@readGridCutter, shape " << name << " with faces not found in " << fileName
              << msg::normal() << '\n';
    return false;
  }
  C.name = name;
  return true;
}

// Number of faces of the cutter crossed by the ray from p along dir (Moller-Trumbore,
// the polygonal faces being split into triangles from their first vertex)
size_t gridCutterCrossings(const GridCutter& C, const vec3r& p, const vec3r& dir) {
  size_t count = 0;
  for (const auto& f : C.face) {
    for (size_t k = 1; k + 1 < f.size(); k++) {
      const vec3r& a = C.vertex[f[0]];
      vec3r e1 = C.vertex[f[k]] - a;
      vec3r e2 = C.vertex[f[k + 1]] - a;
      vec3r h = cross(dir, e2);
      double det = e1 * h;
      if (fabs(det) < 1e-18) continue;
      double inv = 1.0 / det;
      vec3r s = p - a;
      double u = (s * h) * inv;
      if (u < 0.0 || u > 1.0) continue;
      vec3r q = cross(s, e1);
      double v = (dir * q) * inv;
      if (v < 0.0 || u + v > 1.0) continue;
      if ((e2 * q) * inv > 0.0) count++;
    }
  }
  return count;
}

// Is p (in the frame of the cutter, at homothety 1) inside? Three rays vote, so that a
// ray going through an edge or a vertex does not decide alone.
bool gridCutterInside(const GridCutter& C, const vec3r& p) {
  static const vec3r dirs[3] = {vec3r(0.8017, 0.4523, 0.3906), vec3r(-0.3172, 0.8761, -0.3631),
                                vec3r(-0.5517, -0.3291, 0.7663)};
  int votes = 0;
  for (int r = 0; r < 3; r++) {
    if (gridCutterCrossings(C, p, dirs[r]) % 2 == 1) votes++;
  }
  return (votes >= 2);
}

// The same as generatePacking_grid (sameCluster false) or generatePacking_grid_clust (true),
// but only the positions inside the cutter, scaled by cutterHomothety, are kept.
int generatePacking_grid_inside(std::ostream& os, const char* name, vec3r& origBox, vec3r& boxSize, vec3i& n, int group,
                                int clustID, double homothety, int randQ, bool sameCluster, const GridCutter& cutter,
                                double cutterHomothety) {
  int Np = 0;
  quat Qclust;
  vec3r step(boxSize.x / (double)n.x, boxSize.y / (double)n.y, boxSize.z / (double)n.z);
  vec3r orig = origBox + vec3r(0.5 * step.x, 0.5 * step.y, 0.5 * step.z);

  for (int iy = 0; iy < n.y; iy++) {
    for (int ix = 0; ix < n.x; ix++) {
      for (int iz = 0; iz < n.z; iz++) {
        vec3r gridPos = orig + vec3r(ix * step.x, iy * step.y, iz * step.z);
        if (!gridCutterInside(cutter, gridPos / cutterHomothety)) continue;
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
        if (!sameCluster) clustID++;
        Np++;
      }
    }
  }
  return Np;
}

#endif /* end of include guard: GENERATEPACKING_GRID_INSIDE_HPP */
