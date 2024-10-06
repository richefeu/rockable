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
//
// Convert STL to shapeSurvey format
// Author: Vincent.Richefeu@univ-grenoble-alpes.fr
// Lab 3SR, Grenoble
// 2016, 2024

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

#include <map>
#include <set>
#include <vector>

#include <sstream>
#include <stdint.h>
#include <string>

#define ZERO 1e-12

typedef double real;
typedef unsigned int uint;
typedef float real32;
typedef double real64;

#include "tclap/CmdLine.h"

#include "AABB.hpp"
#include "message.hpp"
#include "vec3.hpp"

// ==================
// = DATA STRUCTURE =
// ==================

// point
struct point {
  double x, y, z;  // coordinates
  point() : x(0.0), y(0.0), z(0.0) {}
};

// triangle
struct triangle {
  int i0, i1, i2;  // ID-numbers of points in the mesh
  vec3r normal;
  triangle() : i0(0), i1(0), i2(0), normal() {}
  triangle(int I0, int I1, int I2) : i0(I0), i1(I1), i2(I2), normal() {}
};

// edge
struct edge {
  int i0, i1;  // ID-numbers of points in the mesh
  edge() : i0(0), i1(0) {}
  edge(int I0, int I1) : i0(I0), i1(I1) {}
};

// This specialization of less is useful for sorting edges in lexicographic order
namespace std {
template <>
struct less<edge> {
  bool operator()(const edge& lhs, const edge& rhs) const {
    if (lhs.i0 < rhs.i0) return true;
    if ((lhs.i0 == rhs.i0) && (lhs.i1 < rhs.i1)) return true;
    return false;
  }
};
}  // namespace std

struct mesh {
  std::string name;
  double radius;
  std::vector<point> points;
  std::vector<triangle> triangles;
  std::vector<edge> edges;

  // Ctor
  mesh(const char* Name, double R) : name(Name), radius(R) {}

  void clear() {
    points.clear();
    triangles.clear();
    edges.clear();
  }

  void getAABB(AABB& box) {
    box.min.set(1e20, 1e20, 1e20);
    box.max.set(-1e20, -1e20, -1e20);
    for (size_t i = 0; i < points.size(); i++) {
      if (points[i].x < box.min.x) box.min.x = points[i].x;
      if (points[i].y < box.min.y) box.min.y = points[i].y;
      if (points[i].z < box.min.z) box.min.z = points[i].z;
      if (points[i].x > box.max.x) box.max.x = points[i].x;
      if (points[i].y > box.max.y) box.max.y = points[i].y;
      if (points[i].z > box.max.z) box.max.z = points[i].z;
    }
    box.enlarge(radius);
  }

  void clean() {
    std::set<std::pair<size_t, size_t> > edgeSet;
    for (size_t e = 0; e < edges.size(); e++) {
      std::pair<size_t, size_t> P;
      if (edges[e].i0 < edges[e].i1) {
        P.first = edges[e].i0;
        P.second = edges[e].i1;
      } else {
        P.first = edges[e].i1;
        P.second = edges[e].i0;
      }
      edgeSet.insert(P);  // That way, duplication is not allowed
    }

    size_t nbePrev = edges.size();
    size_t nbeNew = edgeSet.size();

    edges.clear();
    for (auto it : edgeSet) {
      edges.push_back(edge(it.first, it.second));
    }

    std::cout << nbePrev - nbeNew << " duplicated edges have been removed\n";
  }
};

// UINT8[80]  Header
// UINT32     Number of triangles
//
// foreach triangle
// REAL32[3]  Normal vector
// REAL32[3]  Vertex 1
// REAL32[3]  Vertex 2
// REAL32[3]  Vertex 3
// UINT16     Attribute byte count
// end
void readStlBin(const char* name, mesh& Mesh) {
  std::ifstream file(name, std::ifstream::binary | std::ifstream::in);
  if (!file) {
    std::cerr << "Cannot read " << name << std::endl;
    return;
  }

  std::cerr << "Read mesh file (stl BINARY format)... " << std::flush;

  uint8_t head[80];
  file.read((char*)&head, sizeof(uint8_t) * 80);

  uint32_t nb_triangles;
  file.read((char*)&nb_triangles, sizeof(uint32_t));

  vec3<real32> P1, P2, P3;
  triangle T;
  std::map<vec3<real32>, size_t, std::less<vec3<real32> > > map_points;
  std::map<vec3<real32>, size_t, std::less<vec3<real32> > >::iterator itMap;
  size_t ivertex = 0;
  uint16_t attByte;
  int swp;

  Mesh.triangles.clear();
  Mesh.points.clear();

  vec3<real32> normal, vertex1, vertex2, vertex3;

  for (uint32_t t = 0; t < nb_triangles; t++) {
    file.read((char*)&normal, sizeof(real32) * 3);
    file.read((char*)&vertex1, sizeof(real32) * 3);
    file.read((char*)&vertex2, sizeof(real32) * 3);
    file.read((char*)&vertex3, sizeof(real32) * 3);
    file.read((char*)&attByte, sizeof(uint16_t));

    T.normal.x = (real)(normal.x);
    T.normal.y = (real)(normal.y);
    T.normal.z = (real)(normal.z);

    itMap = map_points.find(vertex1);
    if (itMap != map_points.end()) {
      T.i0 = itMap->second;
    } else {
      map_points[vertex1] = ivertex;
      T.i0 = ivertex;
      ivertex++;
    }

    itMap = map_points.find(vertex2);
    if (itMap != map_points.end()) {
      T.i1 = itMap->second;
    } else {
      map_points[vertex2] = ivertex;
      T.i1 = ivertex;
      ivertex++;
    }

    itMap = map_points.find(vertex3);
    if (itMap != map_points.end()) {
      T.i2 = itMap->second;
    } else {
      map_points[vertex3] = ivertex;
      T.i2 = ivertex;
      ivertex++;
    }

    if (T.i0 > T.i1) {
      swp = T.i0;
      T.i0 = T.i1;
      T.i1 = swp;
    }
    if (T.i1 > T.i2) {
      swp = T.i1;
      T.i1 = T.i2;
      T.i2 = swp;
    }
    if (T.i0 > T.i1) {
      swp = T.i0;
      T.i0 = T.i1;
      T.i1 = swp;
    }

    Mesh.triangles.push_back(T);
  }

  std::map<size_t, vec3<real32> > map_id;
  for (itMap = map_points.begin(); itMap != map_points.end(); ++itMap) {
    map_id[itMap->second] = itMap->first;
  }
  map_points.clear();

  point Pt;
  for (std::map<size_t, vec3<real32> >::iterator it = map_id.begin(); it != map_id.end(); ++it) {
    Pt.x = (real)((it->second).x);
    Pt.y = (real)((it->second).y);
    Pt.z = (real)((it->second).z);
    Mesh.points.push_back(Pt);
  }

  // Build the set of edges
  Mesh.edges.clear();
  std::set<edge, std::less<edge> > tmp_edges;
  edge E;
  for (size_t t = 0; t < Mesh.triangles.size(); t++) {
    E.i0 = Mesh.triangles[t].i0;
    E.i1 = Mesh.triangles[t].i1;
    tmp_edges.insert(E);
    E.i0 = Mesh.triangles[t].i1;
    E.i1 = Mesh.triangles[t].i2;
    tmp_edges.insert(E);
    E.i0 = Mesh.triangles[t].i2;
    E.i1 = Mesh.triangles[t].i0;
    tmp_edges.insert(E);
  }

  for (std::set<edge, std::less<edge> >::iterator it = tmp_edges.begin(); it != tmp_edges.end(); ++it) {
    Mesh.edges.push_back(*it);
  }

  std::cout << "done." << std::endl;
  std::cout << "  Number of vertices  " << Mesh.points.size() << std::endl;
  std::cout << "  Number of edges     " << Mesh.edges.size() << std::endl;
  std::cout << "  Number of triangles " << Mesh.triangles.size() << std::endl;
}

// Exportation in shp format (for Rockable)
void exportShape(mesh& Mesh) {
  std::ofstream file("mesh.shp");
  file << "<" << std::endl;
  file << "name " << Mesh.name << std::endl;
  file << "radius " << Mesh.radius << std::endl;
  file << "preCompDone n" << std::endl << std::endl;

  AABB box;
  Mesh.getAABB(box);
  file << "obb.center " << 0.5 * (box.max + box.min) << std::endl;
  file << "obb.extent " << 0.5 * (box.max - box.min) << std::endl;
  file << "obb.e1 1 0 0" << std::endl;
  file << "obb.e2 0 1 0" << std::endl;
  file << "obb.e3 0 0 1" << std::endl;

  file << "nv " << Mesh.points.size() << std::endl;

  for (size_t v = 0; v < Mesh.points.size(); v++) {
    file << Mesh.points[v].x << " " << Mesh.points[v].y << " " << Mesh.points[v].z << std::endl;
  }
  file << std::endl;

  file << "ne " << Mesh.edges.size() << std::endl;
  for (size_t e = 0; e < Mesh.edges.size(); e++) {
    file << Mesh.edges[e].i0 << " " << Mesh.edges[e].i1 << std::endl;
  }
  file << std::endl;

  file << "nf " << Mesh.triangles.size() << std::endl;
  for (size_t f = 0; f < Mesh.triangles.size(); f++) {
    file << "3 " << Mesh.triangles[f].i0 << " " << Mesh.triangles[f].i1 << " " << Mesh.triangles[f].i2 << std::endl;
  }

  file << ">" << std::endl;
}

// This function resize the model so that the largest distance between 2 vertices is 'WantedSize'.
// The radius is not accounted for in the sieving size
void setSieveSize(mesh& Mesh, double WantedSize, bool scaleRadius = true) {
  double MaxSize = 0.0;

  double s = 0.0;
  vec3r l;
  for (size_t i = 0; i < Mesh.points.size(); i++) {
    for (size_t j = i + 1; j < Mesh.points.size(); j++) {
      l.set(Mesh.points[j].x - Mesh.points[i].x, Mesh.points[j].y - Mesh.points[i].y,
            Mesh.points[j].z - Mesh.points[i].z);
      s = l.length();
      if (s > MaxSize) MaxSize = s;
    }
  }

  if (MaxSize == 0.0) return;

  double scaleFactor = WantedSize / MaxSize;
  for (size_t i = 0; i < Mesh.points.size(); i++) {
    Mesh.points[i].x *= scaleFactor;
    Mesh.points[i].y *= scaleFactor;
    Mesh.points[i].z *= scaleFactor;
  }
  if (scaleRadius == true) {
    Mesh.radius *= scaleFactor;
    std::cout << "radius has been rescaled to " << Mesh.radius << '\n';
  }
}

// This function rescale the model (and possibly the radius)
void rescale(mesh& Mesh, double scaleFactor, bool scaleRadius = true) {
  for (size_t i = 0; i < Mesh.points.size(); i++) {
    Mesh.points[i].x *= scaleFactor;
    Mesh.points[i].y *= scaleFactor;
    Mesh.points[i].z *= scaleFactor;
  }
  if (scaleRadius == true) {
    Mesh.radius *= scaleFactor;
    std::cout << "radius has been rescaled to " << Mesh.radius << '\n';
  }
}

int main(int argc, char const* argv[]) {
  std::string stlFileName;
  double scaleFactor = 1.0;
  double maxLength = 0.0;
  double radius = 1.0;
  bool clean = false;
  bool scaleRadius = false;

  try {
    TCLAP::CmdLine cmd("Convert a binary STL file to a shape that can be used by Rockable", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the input STL binary file", true, "file.stl",
                                                  "stl file");
    TCLAP::ValueArg<double> radiusArg("r", "radius", "Radius of the rounded edges (Minskowski radius) to be used", true,
                                      1.0, "double");
    TCLAP::ValueArg<double> scaleArg("s", "scaleFactor", "Set the (re-)scale factor of the output object", false, 1.0,
                                     "double");
    TCLAP::ValueArg<double> maxLengthArg("m", "maxLength", "Set the maximum length (sieving size) of the output object",
                                         false, 0.0, "double");
    TCLAP::SwitchArg cleanArg("c", "clean", "Remove duplicated edges", false);
    TCLAP::SwitchArg scaleRadiusArg("z", "scaleRadius", "Rescale the radius", false);

    cmd.add(nameArg);
    cmd.add(radiusArg);
    cmd.add(scaleArg);
    cmd.add(maxLengthArg);
    cmd.add(cleanArg);
    cmd.add(scaleRadiusArg);

    cmd.parse(argc, argv);

    stlFileName = nameArg.getValue();
    scaleFactor = scaleArg.getValue();
    maxLength = maxLengthArg.getValue();
    radius = radiusArg.getValue();
    clean = cleanArg.getValue();
    scaleRadius = scaleRadiusArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  mesh MESH(stlFileName.c_str(), radius);
  readStlBin(stlFileName.c_str(), MESH);

  if (clean == true) {
    std::cout << "Remove duplicated edges\n";
    MESH.clean();
    std::cout << "  Number of edges     " << MESH.edges.size() << '\n';
  }

  if (maxLength > 0.0) {
    std::cout << "Set sieve size to " << maxLength << '\n';
    setSieveSize(MESH, maxLength, scaleRadius);
  }

  if (scaleFactor != 1.0) {
    std::cout << "Rescale by factor " << scaleFactor << '\n';
    rescale(MESH, scaleFactor, scaleRadius);
  }

  exportShape(MESH);

  return 0;
}
