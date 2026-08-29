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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <unordered_map>

#include "AABB.hpp"
#include "ShapeMesh.hpp"
#include "marchingCubesTables.hpp"

double ShapeMesh::computeVolume() const {
  // Signed volume of a closed triangle mesh: sum of the signed volumes of the
  // tetrahedra (origin, a, b, c). Robust regardless of where the origin sits
  double vol = 0.0;
  for (size_t t = 0; t < tri.size(); ++t) {
    const vec3r& a = P[tri[t].a];
    const vec3r& b = P[tri[t].b];
    const vec3r& c = P[tri[t].c];
    vol += (a * cross(b, c)) / 6.0;
  }
  return std::fabs(vol);
}

double ShapeMesh::computeArea() const {
  double area = 0.0;
  for (size_t t = 0; t < tri.size(); ++t) {
    const vec3r& a = P[tri[t].a];
    const vec3r& b = P[tri[t].b];
    const vec3r& c = P[tri[t].c];
    area += 0.5 * norm(cross(b - a, c - a));
  }
  return area;
}

void ShapeMesh::writeOBJ(std::ostream& os) const {
  os << "# r-shape skin mesh (Rockable)\n";
  os << "# vertices " << P.size() << "  triangles " << tri.size() << "  volume " << volume << "  area " << area << "\n";
  for (size_t i = 0; i < P.size(); ++i) {
    os << "v " << P[i].x << ' ' << P[i].y << ' ' << P[i].z << '\n';
  }
  for (size_t i = 0; i < N.size(); ++i) {
    os << "vn " << N[i].x << ' ' << N[i].y << ' ' << N[i].z << '\n';
  }
  // OBJ indices are 1-based; positions and normals share the same index here
  for (size_t t = 0; t < tri.size(); ++t) {
    size_t a = tri[t].a + 1, b = tri[t].b + 1, c = tri[t].c + 1;
    os << "f " << a << "//" << a << ' ' << b << "//" << b << ' ' << c << "//" << c << '\n';
  }
}

void ShapeMesh::writePLY(std::ostream& os) const {
  os << "ply\n";
  os << "format ascii 1.0\n";
  os << "comment r-shape skin mesh (Rockable), volume " << volume << " area " << area << "\n";
  os << "element vertex " << P.size() << "\n";
  os << "property float x\nproperty float y\nproperty float z\n";
  os << "property float nx\nproperty float ny\nproperty float nz\n";
  os << "element face " << tri.size() << "\n";
  os << "property list uchar int vertex_indices\n";
  os << "end_header\n";
  for (size_t i = 0; i < P.size(); ++i) {
    os << P[i].x << ' ' << P[i].y << ' ' << P[i].z << ' ' << N[i].x << ' ' << N[i].y << ' ' << N[i].z << '\n';
  }
  for (size_t t = 0; t < tri.size(); ++t) {
    os << "3 " << tri[t].a << ' ' << tri[t].b << ' ' << tri[t].c << '\n';
  }
}

// One block of the ".rmsh" companion file: a name, the vertices (position and
// exact normal) and the triangles. Kept simple and text-based, in the spirit of
// the ".shp" shape file it sits next to.
void ShapeMesh::writeRmsh(std::ostream& os, const std::string& name) const {
  os << "<\n";
  os << "name " << name << "\n";
  os << "nv " << P.size() << "\n";
  for (size_t i = 0; i < P.size(); ++i) {
    os << P[i].x << ' ' << P[i].y << ' ' << P[i].z << ' ' << N[i].x << ' ' << N[i].y << ' ' << N[i].z << '\n';
  }
  os << "nt " << tri.size() << "\n";
  for (size_t t = 0; t < tri.size(); ++t) {
    os << tri[t].a << ' ' << tri[t].b << ' ' << tri[t].c << '\n';
  }
  os << ">\n";
}

std::map<std::string, ShapeMesh> loadRmsh(const std::string& path) {
  std::map<std::string, ShapeMesh> meshes;
  std::ifstream is(path.c_str());
  if (!is) return meshes;

  std::string tok;
  while (is >> tok) {
    if (tok != "<") continue;  // start of a shape block
    ShapeMesh m;
    std::string name;
    while (is >> tok && tok != ">") {
      if (tok == "name") {
        is >> name;
      } else if (tok == "nv") {
        size_t nv = 0;
        is >> nv;
        m.P.resize(nv);
        m.N.resize(nv);
        for (size_t i = 0; i < nv; ++i) {
          is >> m.P[i].x >> m.P[i].y >> m.P[i].z >> m.N[i].x >> m.N[i].y >> m.N[i].z;
        }
      } else if (tok == "nt") {
        size_t nt = 0;
        is >> nt;
        m.tri.resize(nt);
        for (size_t t = 0; t < nt; ++t) {
          is >> m.tri[t].a >> m.tri[t].b >> m.tri[t].c;
        }
      }
    }
    if (!name.empty()) meshes[name] = m;
  }
  return meshes;
}

namespace {

// The eight corners of a marching-cubes cell, in the order the classic tables
// expect (Paul Bourke's convention, see marchingCubesTables.hpp)
const int CORNER[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                          {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// The twelve edges of a cell, as pairs of corner indices, same convention
const int EDGE_CORNER[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

struct GridSampler {
  const ShapeSDF* sdf;
  vec3r origin;
  double h;
  int nx, ny, nz;  // number of grid nodes along each axis
  std::vector<double> val;

  double& at(int i, int j, int k) { return val[(size_t)(k * ny + j) * nx + i]; }

  vec3r nodePos(int i, int j, int k) const {
    return origin + vec3r(i * h, j * h, k * h);
  }

  void sampleAll() {
    val.assign((size_t)nx * ny * nz, 0.0);
    // Every node is independent. This is the dominant cost on a detailed shape,
    // and it parallelises perfectly (the pragma is simply ignored without OpenMP)
#pragma omp parallel for schedule(static)
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          at(i, j, k) = sdf->value(nodePos(i, j, k));
        }
      }
    }
  }
};

// A grid edge is identified by its lower corner (i,j,k) and its axis (0,1,2).
// This gives every edge one global key, so that a marching-cubes vertex placed
// on a shared edge is created once and the mesh comes out watertight
uint64_t edgeKey(int i, int j, int k, int axis) {
  uint64_t key = (uint64_t)(uint32_t)i;
  key = key * 2048 + (uint64_t)(uint32_t)j;
  key = key * 2048 + (uint64_t)(uint32_t)k;
  key = key * 4 + (uint64_t)axis;
  return key;
}

}  // namespace

ShapeMesh buildShapeMesh(const Shape& shape, const ShapeMeshOptions& opt) {
  ShapeMesh mesh;
  ShapeSDF sdf(shape);

  double R = shape.radius;

  // Target sag epsilon, then the grid step from the sphere sag relation
  // sag ~ R * theta^2 / 8 with chord h ~ R * theta, hence h ~ sqrt(8 eps R).
  // A thin shell has thickness 2R, so cap h at R/2 to never step over it
  double epsilon = (opt.epsilon > 0.0) ? opt.epsilon : opt.epsilonRel * R;
  if (epsilon <= 0.0) epsilon = 1.0e-2;  // last resort, e.g. radius == 0
  mesh.epsilon = epsilon;

  double h = std::sqrt(8.0 * epsilon * std::max(R, epsilon));
  // A shape that is not a thick solid (a sphere, capsule, plate or open surface)
  // is entirely at the R scale, or a slab of thickness 2R; the grid must not
  // step over it and make it vanish. A filled solid has a thick interior, so
  // capping its step here would only make it coarse and epsilon-insensitive
  if (!sdf.isSolid() && R > 0.0) h = std::min(h, 0.5 * R);
  if (h <= 0.0) h = epsilon;

  AABB box;
  const_cast<Shape&>(shape).getAABB(box);  // getAABB already enlarges by radius
  vec3r span = box.max - box.min;

  // Hard budget on the number of grid nodes, so a large thin shape (a wall
  // plate, say) cannot blow up into a runaway grid. If the requested step would
  // exceed it, coarsen h to fit and warn: the mesh is then cruder than asked
  if (opt.maxCells > 0) {
    double n0 = (std::ceil(span.x / h) + 1 + 2 * opt.marginCells);
    double n1 = (std::ceil(span.y / h) + 1 + 2 * opt.marginCells);
    double n2 = (std::ceil(span.z / h) + 1 + 2 * opt.marginCells);
    if (n0 * n1 * n2 > (double)opt.maxCells) {
      double factor = std::cbrt((n0 * n1 * n2) / (double)opt.maxCells);
      double hClamped = h * factor;
      std::fprintf(stderr,
                   "[ShapeMesh] '%s': grid capped at %.0f nodes; step coarsened %.2fx "
                   "(sag ~ %.3g instead of %.3g). Pass a larger epsilon or maxCells to refine.\n",
                   shape.name.c_str(), (double)opt.maxCells, factor, hClamped * hClamped / (8.0 * std::max(R, epsilon)),
                   epsilon);
      h = hClamped;
    }
  }

  GridSampler g;
  g.sdf = &sdf;
  g.h = h;
  g.origin = box.min - vec3r(opt.marginCells * h, opt.marginCells * h, opt.marginCells * h);
  g.nx = (int)std::ceil(span.x / h) + 1 + 2 * opt.marginCells;
  g.ny = (int)std::ceil(span.y / h) + 1 + 2 * opt.marginCells;
  g.nz = (int)std::ceil(span.z / h) + 1 + 2 * opt.marginCells;
  g.sampleAll();

  std::unordered_map<uint64_t, size_t> edgeVertex;

  // Marching cubes over every cell
  for (int k = 0; k + 1 < g.nz; ++k) {
    for (int j = 0; j + 1 < g.ny; ++j) {
      for (int i = 0; i + 1 < g.nx; ++i) {
        double cval[8];
        int cubeIndex = 0;
        for (int c = 0; c < 8; ++c) {
          cval[c] = g.at(i + CORNER[c][0], j + CORNER[c][1], k + CORNER[c][2]);
          if (cval[c] < 0.0) cubeIndex |= (1 << c);
        }

        int edges = MC_EDGE_TABLE[cubeIndex];
        if (edges == 0) continue;

        size_t edgeVtx[12];
        for (int e = 0; e < 12; ++e) {
          if (!(edges & (1 << e))) continue;

          int c0 = EDGE_CORNER[e][0];
          int c1 = EDGE_CORNER[e][1];
          int i0 = i + CORNER[c0][0], j0 = j + CORNER[c0][1], k0 = k + CORNER[c0][2];
          int i1 = i + CORNER[c1][0], j1 = j + CORNER[c1][1], k1 = k + CORNER[c1][2];

          // Global key of this grid edge, taken from its lower endpoint and axis
          int axis = (i1 != i0) ? 0 : ((j1 != j0) ? 1 : 2);
          int li = std::min(i0, i1), lj = std::min(j0, j1), lk = std::min(k0, k1);
          uint64_t key = edgeKey(li, lj, lk, axis);

          std::unordered_map<uint64_t, size_t>::iterator it = edgeVertex.find(key);
          if (it != edgeVertex.end()) {
            edgeVtx[e] = it->second;
            continue;
          }

          // Linear guess along the edge from the two corner values, then refine
          double f0 = cval[c0];
          double f1 = cval[c1];
          double t = (std::fabs(f1 - f0) > 1.0e-30) ? (f0 / (f0 - f1)) : 0.5;
          vec3r p = g.nodePos(i0, j0, k0) + t * (g.nodePos(i1, j1, k1) - g.nodePos(i0, j0, k0));

          // Project onto f = 0 by Newton along the analytic gradient
          for (int nit = 0; nit < opt.newtonIters; ++nit) {
            double f = sdf.value(p);
            vec3r grad = sdf.gradient(p);
            if (grad.isnull()) break;
            p -= f * grad;  // grad is unit length, so this is a full Newton step
          }

          size_t idx = mesh.P.size();
          mesh.P.push_back(p);
          vec3r n = sdf.gradient(p);
          mesh.N.push_back(n);
          ClosestPoint cp = sdf.closestPoint(p);
          mesh.primType.push_back(cp.primType);
          mesh.primIndex.push_back(cp.primIndex);

          edgeVertex[key] = idx;
          edgeVtx[e] = idx;
        }

        const int8_t* row = MC_TRI_TABLE[cubeIndex];
        for (int t = 0; row[t] != -1; t += 3) {
          Tri tr;
          tr.a = edgeVtx[row[t]];
          tr.b = edgeVtx[row[t + 1]];
          tr.c = edgeVtx[row[t + 2]];
          mesh.tri.push_back(tr);
        }
      }
    }
  }

  mesh.volume = mesh.computeVolume();
  mesh.area = mesh.computeArea();
  return mesh;
}
