// Shape builders shared by the ShapeSDF and ShapeMesh unit tests.

#ifndef TEST_SHAPES_HPP
#define TEST_SHAPES_HPP

#include <fstream>
#include <string>
#include <vector>

#include "Core/Shape.hpp"

// A single vertex: a sphere of radius R
inline Shape makeSphere(double R) {
  Shape s;
  s.name = "sphere";
  s.radius = R;
  s.vertex.push_back(vec3r(0, 0, 0));
  return s;
}

// Two vertices and one edge: a capsule of radius R and axis length L
inline Shape makeCapsule(double R, double L) {
  Shape s;
  s.name = "capsule";
  s.radius = R;
  s.vertex.push_back(vec3r(-0.5 * L, 0, 0));
  s.vertex.push_back(vec3r(0.5 * L, 0, 0));
  s.edge.push_back(std::pair<size_t, size_t>(0, 1));
  return s;
}

// A square plate of half-side a in z = 0, declared as an open surface
inline Shape makePlate(double R, double a) {
  Shape s;
  s.name = "plate";
  s.radius = R;
  s.isSurface = true;
  s.vertex.push_back(vec3r(-a, -a, 0));
  s.vertex.push_back(vec3r(a, -a, 0));
  s.vertex.push_back(vec3r(a, a, 0));
  s.vertex.push_back(vec3r(-a, a, 0));
  for (size_t i = 0; i < 4; ++i) s.edge.push_back(std::pair<size_t, size_t>(i, (i + 1) % 4));
  std::vector<size_t> f;
  for (size_t i = 0; i < 4; ++i) f.push_back(i);
  s.face.push_back(f);
  return s;
}

// Cube of half-side a. Faces are fed with a deliberately inconsistent winding,
// to exercise the outward re-orientation
inline Shape makeCube(double R, double a) {
  Shape s;
  s.name = "cube";
  s.radius = R;
  const double v[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                          {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  for (size_t i = 0; i < 8; ++i) s.vertex.push_back(vec3r(a * v[i][0], a * v[i][1], a * v[i][2]));
  const size_t e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                           {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  for (size_t i = 0; i < 12; ++i) s.edge.push_back(std::pair<size_t, size_t>(e[i][0], e[i][1]));
  const size_t f[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  for (size_t i = 0; i < 6; ++i) s.face.push_back(std::vector<size_t>(f[i], f[i] + 4));
  return s;
}

// A cube that is NOT a clean 2-manifold: same solid as makeCube, plus one
// degenerate (2-vertex) face. That extra face makes the manifold test fail, so
// the SDF must fall back to ray casting for the sign, exactly as it does on an
// imperfect STL soup, while still meshing the correct solid
inline Shape makeCubeSoup(double R, double a) {
  Shape s = makeCube(R, a);
  s.name = "cubeSoup";
  std::vector<size_t> bogus;  // a face with only two vertices: geometrically void
  bogus.push_back(0);
  bogus.push_back(6);
  s.face.push_back(bogus);
  return s;
}

// A concave L-shaped prism: an L polygon extruded along z. The vertical edge at
// the notch is a reentrant (concave) edge, where the dilation keeps a sharp
// crease and where the sign relies on the angle-weighted pseudonormal. The L
// polygon (area 3 = a 2x2 square minus a 1x1 corner) is centered on the origin
inline Shape makeL(double R, double height = 1.0) {
  Shape s;
  s.name = "L";
  s.radius = R;
  // L polygon, counter-clockwise, with the reentrant corner at (1,1)
  const double p[6][2] = {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}};
  double cx = 1.0, cy = 1.0;  // recenter roughly
  double zb = -0.5 * height, zt = 0.5 * height;
  for (int i = 0; i < 6; ++i) s.vertex.push_back(vec3r(p[i][0] - cx, p[i][1] - cy, zb));
  for (int i = 0; i < 6; ++i) s.vertex.push_back(vec3r(p[i][0] - cx, p[i][1] - cy, zt));

  for (int i = 0; i < 6; ++i) {  // bottom ring, top ring, verticals
    s.edge.push_back(std::pair<size_t, size_t>(i, (i + 1) % 6));
    s.edge.push_back(std::pair<size_t, size_t>(6 + i, 6 + (i + 1) % 6));
    s.edge.push_back(std::pair<size_t, size_t>(i, 6 + i));
  }

  std::vector<size_t> bottom, top;
  for (int i = 0; i < 6; ++i) bottom.push_back(i);
  for (int i = 0; i < 6; ++i) top.push_back(6 + i);
  s.face.push_back(bottom);
  s.face.push_back(top);
  for (int i = 0; i < 6; ++i) {  // one quad per bottom edge
    int j = (i + 1) % 6;
    std::vector<size_t> side;
    side.push_back(i);
    side.push_back(j);
    side.push_back(6 + j);
    side.push_back(6 + i);
    s.face.push_back(side);
  }
  return s;
}

// Read the first shape ('<' ... '>' block) found in a .shp file. Returns an
// empty (nameless) shape if the file cannot be opened
inline Shape loadFirstShape(const std::string& path) {
  Shape s;
  std::ifstream is(path.c_str());
  if (!is) return s;
  std::string tok;
  while (is >> tok) {
    if (tok == "<") {  // Shape::read consumes the keywords up to the closing '>'
      s.read(is);
      break;
    }
  }
  return s;
}

#endif /* end of include guard: TEST_SHAPES_HPP */
