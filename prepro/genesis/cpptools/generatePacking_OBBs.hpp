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

#ifndef GENERATEPACKING_OBBS_HPP
#define GENERATEPACKING_OBBS_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "OBB.hpp"
#include "message.hpp"
#include "quat.hpp"
#include "transformation.hpp"
#include "vec3.hpp"

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp

// Random sequential packing of oriented bounding boxes (OBB) in a container,
// in the way of TubePacker: the boxes are placed one by one, the largest first,
// at a random position and with a random orientation, until they no longer
// overlap the boxes already placed. Each box is then replaced by the particle
// it bounds, so that a computation can start from a dense state, without any
// overlap between the particles (it still has to be stabilised).
//
// The OBB of a shape is read in the shape file (obb.extent, obb.center, obb.e1,
// obb.e2, obb.e3): it includes the Minkowski radius, and it is not necessarily
// centred on the mass centre of the shape, nor aligned with its principal axes.

struct OBBPackShape {
  std::string name;
  vec3r center;        // centre of the OBB in the body frame
  vec3r e[3];          // axes of the OBB in the body frame
  vec3r extent;        // half-sizes of the OBB
  double volume{0.0};  // volume of the shape (0 if unknown)
};

// Reads the OBB of the shape 'name' in the shape file 'fileName'.
// Returns false (with a message) if the shape is not found or not pre-computed.
bool readOBBPackShape(const char* fileName, const std::string& name, OBBPackShape& S) {
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::cerr << msg::warn() << "@readOBBPackShape, cannot read " << fileName << msg::normal() << '\n';
    return false;
  }

  std::string line;
  bool inShape = false;
  bool found = false;
  bool preCompDone = false;
  bool hasOBB = false;
  while (std::getline(file, line)) {
    std::istringstream is(line);
    std::string key;
    if (!(is >> key)) continue;
    if (key == "<") {
      inShape = false;
      preCompDone = false;
      hasOBB = false;
      S = OBBPackShape();
      S.e[0].set(1.0, 0.0, 0.0);
      S.e[1].set(0.0, 1.0, 0.0);
      S.e[2].set(0.0, 0.0, 1.0);
    } else if (key == ">") {
      if (inShape) {
        found = true;
        break;
      }
    } else if (key == "name") {
      std::string n;
      is >> n;
      inShape = (n == name);
    } else if (key == "nv" || key == "ne" || key == "nf") {
      size_t n = 0;
      is >> n;
      for (size_t i = 0; i < n; i++) std::getline(file, line);
    } else if (key == "preCompDone") {
      std::string yn;
      is >> yn;
      preCompDone = (yn == "y" || yn == "Y");
    } else if (key == "obb.extent") {
      is >> S.extent;
      hasOBB = true;
    } else if (key == "obb.center") {
      is >> S.center;
    } else if (key == "obb.e1") {
      is >> S.e[0];
    } else if (key == "obb.e2") {
      is >> S.e[1];
    } else if (key == "obb.e3") {
      is >> S.e[2];
    } else if (key == "volume") {
      is >> S.volume;
    }
  }

  if (!found) {
    std::cerr << msg::warn() << "@readOBBPackShape, shape " << name << " not found in " << fileName << msg::normal()
              << '\n';
    return false;
  }
  if (!preCompDone || !hasOBB) {
    std::cerr << msg::warn() << "@readOBBPackShape, the OBB of shape " << name << " is not computed: run 'shapeSurvey "
              << fileName << " -P " << fileName << "' first" << msg::normal() << '\n';
    return false;
  }
  S.name = name;
  for (int k = 0; k < 3; k++) S.e[k].normalize();
  return true;
}

struct OBBPackItem {
  size_t shape{0};
  double homothety{1.0};
  double boundingRadius{0.0};
  double obbVolume{0.0};
  quat Q;   // orientation of the body
  OBB obb;  // the OBB in the global frame
};

// Uniformly distributed random orientation (Shoemake)
template <typename Engine>
quat OBBPack_randomQuat(Engine& engine) {
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  double u1 = u01(engine);
  double u2 = 2.0 * M_PI * u01(engine);
  double u3 = 2.0 * M_PI * u01(engine);
  double a = sqrt(1.0 - u1);
  double b = sqrt(u1);
  return quat(a * sin(u2), a * cos(u2), b * sin(u3), b * cos(u3));
}

// The ways of placing the boxes
const int OBBPACK_RANDOM = 0;  // anywhere in the container (random sequential addition)
const int OBBPACK_STACK = 1;   // as low as possible, along the axis of the container

// Packs OBBs in a container, and writes the particles they bound.
//
// - containerType: 0 = cuboid, 1 = cylinder (elliptic if its two transverse sizes differ)
// - mode: OBBPACK_RANDOM, each box is put at a random place and with a random orientation, until it
//   overlaps no other box (random sequential addition, as TubePacker); OBBPACK_STACK, each box is
//   tried 'nStackTrials' times, with a random orientation and a random transverse position, and
//   goes down along the axis to the lowest place where it overlaps no other box; the lowest trial
//   is kept. Stacking gives denser packings, with the particles at the bottom of the container.
// - axis: 0, 1 or 2 (x, y or z), the axis of the cylinder, the direction of stacking, and the
//   direction along which the particles are sorted (the first ones written are at the bottom)
// - origBox, boxSize: the corner and the sizes of the container (the bounding box of the cylinder)
//   If boxSize along 'axis' is 0, the container has no top: with OBBPACK_RANDOM, its height is
//   computed so that the nb boxes reach the solid fraction 'solidFraction'; with OBBPACK_STACK,
//   the nb boxes are stacked as high as needed. Otherwise, the packing stops when nb boxes are
//   placed, when the solid fraction of the boxes reaches 'solidFraction', or when a box finds no
//   place.
// - hmin, hmax: the homothety of each particle is drawn uniformly in [hmin, hmax]
// - each particle is given its own cluster, numbered from clustID upwards (if sameCluster is false)
//
// Returns the number of particles written.
int generatePacking_OBBs(std::ostream& os, const std::vector<OBBPackShape>& shapes, int containerType, int mode,
                         int axis, vec3r origBox, vec3r boxSize, double hmin, double hmax, size_t nb,
                         double solidFraction, int group, int clustID, unsigned int seed, bool sameCluster = false,
                         size_t maxTrials = 100000, size_t nStackTrials = 20) {
  if (shapes.empty() || nb == 0) return 0;
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);

  int t1 = (axis + 1) % 3;  // the transverse directions
  int t2 = (axis + 2) % 3;
  double section = boxSize[t1] * boxSize[t2];
  if (containerType == 1) section *= 0.25 * M_PI;

  // Draw the shapes and the homotheties
  std::vector<OBBPackItem> items;
  double sumVobb = 0.0;
  bool noTop = (boxSize[axis] <= 0.0);
  double Vbox = section * boxSize[axis];
  for (size_t i = 0; i < nb; i++) {
    OBBPackItem it;
    it.shape = (shapes.size() == 1) ? 0 : std::min(shapes.size() - 1, (size_t)(u01(engine) * shapes.size()));
    it.homothety = hmin + (hmax - hmin) * u01(engine);
    const OBBPackShape& S = shapes[it.shape];
    double h = it.homothety;
    it.obbVolume = 8.0 * S.extent.x * S.extent.y * S.extent.z * h * h * h;
    it.boundingRadius = h * norm(S.extent);
    items.push_back(it);
    sumVobb += it.obbVolume;
    if (!noTop && sumVobb / Vbox >= solidFraction) break;
  }
  if (noTop && mode == OBBPACK_RANDOM) {
    boxSize[axis] = sumVobb / (solidFraction * section);
    Vbox = section * boxSize[axis];
    noTop = false;
    std::cerr << "@generatePacking:OBBs, height of the container = " << boxSize[axis] << '\n';
  }
  if (mode == OBBPACK_RANDOM) {  // the largest first
    std::sort(items.begin(), items.end(),
              [](const OBBPackItem& a, const OBBPackItem& b) { return a.obbVolume > b.obbVolume; });
  }

  // Cells of the size of the largest bounding sphere, to find the neighbours
  double cellSize = 0.0;
  for (const auto& it : items) cellSize = std::max(cellSize, 2.0 * it.boundingRadius);
  auto cellKey = [&](int cx, int cy, int cz) -> long long {
    return ((long long)(cx + 1000000) * 2000003LL + (long long)(cy + 1000000)) * 2000003LL + (long long)(cz + 1000000);
  };
  auto cellOf = [&](const vec3r& p, int c[3]) {
    for (int d = 0; d < 3; d++) c[d] = (int)floor((p[d] - origBox[d]) / cellSize);
  };
  std::unordered_map<long long, std::vector<size_t>> cells;
  std::vector<OBBPackItem> placed;

  auto overlaps = [&](const OBBPackItem& it) {
    int c[3];
    cellOf(it.obb.center, c);
    for (int ix = c[0] - 1; ix <= c[0] + 1; ix++) {
      for (int iy = c[1] - 1; iy <= c[1] + 1; iy++) {
        for (int iz = c[2] - 1; iz <= c[2] + 1; iz++) {
          auto found = cells.find(cellKey(ix, iy, iz));
          if (found == cells.end()) continue;
          for (size_t j : found->second) {
            double sumR = it.boundingRadius + placed[j].boundingRadius;
            if (norm2(it.obb.center - placed[j].obb.center) > sumR * sumR) continue;
            if (it.obb.intersect(placed[j].obb)) return true;
          }
        }
      }
    }
    return false;
  };

  // Random orientation; returns false if the box is wider than the container
  auto orient = [&](OBBPackItem& it, vec3r& w) {
    const OBBPackShape& S = shapes[it.shape];
    it.Q = OBBPack_randomQuat(engine);
    for (int k = 0; k < 3; k++) it.obb.e[k] = it.Q * S.e[k];
    it.obb.extent = it.homothety * S.extent;
    for (int d = 0; d < 3; d++) {  // half-widths of the box along x, y and z
      w[d] = 0.0;
      for (int k = 0; k < 3; k++) w[d] += it.obb.extent[k] * fabs(it.obb.e[k][d]);
    }
    for (int d = 0; d < 3; d++) {
      if (d == axis && boxSize[d] <= 0.0) continue;
      if (2.0 * w[d] > boxSize[d]) return false;
    }
    return true;
  };

  // The eight corners must be inside the section of the cylinder
  auto insideCylinder = [&](const OBBPackItem& it) {
    if (containerType != 1) return true;
    double a = 0.5 * boxSize[t1];
    double b = 0.5 * boxSize[t2];
    double c1 = origBox[t1] + a;
    double c2 = origBox[t2] + b;
    for (int corner = 0; corner < 8; corner++) {
      vec3r p = it.obb.center;
      for (int k = 0; k < 3; k++) {
        double sign = ((corner >> k) & 1) ? 1.0 : -1.0;
        p += sign * it.obb.extent[k] * it.obb.e[k];
      }
      double x = (p[t1] - c1) / a;
      double y = (p[t2] - c2) / b;
      if (x * x + y * y > 1.0) return false;
    }
    return true;
  };

  double sumVplaced = 0.0;
  double sumVparticles = 0.0;
  double top = origBox[axis];  // the highest point of the boxes
  bool allVolumes = true;
  for (size_t i = 0; i < items.size(); i++) {
    OBBPackItem& it = items[i];
    bool ok = false;

    if (mode == OBBPACK_STACK) {
      OBBPackItem best;
      double bestLevel = 0.0;
      size_t nTried = 0;
      for (size_t trial = 0; trial < maxTrials && nTried < nStackTrials; trial++) {
        vec3r w;
        if (!orient(it, w)) continue;
        for (int d : {t1, t2}) it.obb.center[d] = origBox[d] + w[d] + (boxSize[d] - 2.0 * w[d]) * u01(engine);
        if (!insideCylinder(it)) continue;
        nTried++;
        // go up from the bottom until the box overlaps no other one
        double step = 0.25 * std::min(it.obb.extent.x, std::min(it.obb.extent.y, it.obb.extent.z));
        double level = origBox[axis] + w[axis];
        double levelMax = noTop ? top + w[axis] + step : origBox[axis] + boxSize[axis] - w[axis];
        bool free = false;
        for (; level <= levelMax; level += step) {
          it.obb.center[axis] = level;
          if (!overlaps(it)) {
            free = true;
            break;
          }
        }
        if (free && (!ok || level < bestLevel)) {
          ok = true;
          bestLevel = level;
          best = it;
        }
      }
      if (ok) it = best;
    } else {
      for (size_t trial = 0; trial < maxTrials && !ok; trial++) {
        vec3r w;
        if (!orient(it, w)) continue;
        for (int d = 0; d < 3; d++) it.obb.center[d] = origBox[d] + w[d] + (boxSize[d] - 2.0 * w[d]) * u01(engine);
        if (!insideCylinder(it)) continue;
        if (!overlaps(it)) ok = true;
      }
    }

    if (!ok) {
      std::cerr << "@generatePacking:OBBs, no place found for box " << i + 1 << ", the container is full\n";
      break;
    }

    int c[3];
    cellOf(it.obb.center, c);
    cells[cellKey(c[0], c[1], c[2])].push_back(placed.size());
    placed.push_back(it);
    double halfHeight = 0.0;
    for (int k = 0; k < 3; k++) halfHeight += it.obb.extent[k] * fabs(it.obb.e[k][axis]);
    top = std::max(top, it.obb.center[axis] + halfHeight);
    sumVplaced += it.obbVolume;
    double h = it.homothety;
    if (shapes[it.shape].volume > 0.0) {
      sumVparticles += shapes[it.shape].volume * h * h * h;
    } else {
      allVolumes = false;
    }
  }

  // The particles at the bottom first
  std::sort(placed.begin(), placed.end(),
            [axis](const OBBPackItem& a, const OBBPackItem& b) { return a.obb.center[axis] < b.obb.center[axis]; });

  msg::bestPrecision(os);
  for (size_t i = 0; i < placed.size(); i++) {
    const OBBPackItem& it = placed[i];
    const OBBPackShape& S = shapes[it.shape];
    // The OBB centre is at pos + Q * (h * centre in the body frame)
    vec3r pos = it.obb.center - it.Q * (it.homothety * S.center);
    globalTransformation.apply(pos);
    quat Q = individualParticleRotation * it.Q;
    int clust = sameCluster ? clustID : clustID + (int)i;
    os << S.name << ' ' << group << ' ' << clust << ' ' << it.homothety << "  " << pos << "  0 0 0  0 0 0  " << Q
       << "  0 0 0  0 0 0\n";
  }
  msg::normalPrecision(os);

  // The solid fraction is measured up to the top of the boxes
  double height = noTop ? top - origBox[axis] : boxSize[axis];
  if (mode == OBBPACK_STACK) height = top - origBox[axis];
  double V = section * height;
  std::cerr << "@generatePacking:OBBs, " << placed.size() << " particles, height " << height
            << ", solid fraction of the boxes = " << sumVplaced / V;
  if (allVolumes) std::cerr << ", of the particles = " << sumVparticles / V;
  std::cerr << '\n';

  return (int)placed.size();
}

#endif /* end of include guard: GENERATEPACKING_OBBS_HPP */
