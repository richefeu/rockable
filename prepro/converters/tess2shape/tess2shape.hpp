//
//  tessToshp.hpp
//  rockable
//
//  Created by amarsid on 25/09/2019.
//  Copyright © 2019 amarsid. All rights reserved.
//

#ifndef TESS2SHAPE_HPP
#define TESS2SHAPE_HPP

#include <cstdlib>
#include <vector>

#include "kwParser.hpp"
#include "message.hpp"
#include "polyhTool.hpp"

#include "Core/Shape.hpp"

template <class myType>
int getIndex(std::vector<myType> v, myType b) {
  auto it = find(v.begin(), v.end(), b);
  if (it != v.end()) {
    return int(it - v.begin());
  }
  return -1;
}

struct plan_ {
  vec3r pos;
  vec3r normal;
  double d;
};

struct vertex_ {
  int id{0};
  vec3r v;
  bool operator==(const vertex_& a) { return std::abs(id) == std::abs(a.id); }
  bool operator<(const vertex_& a) { return std::abs(id) < std::abs(a.id); }
};

struct edge_ {
  int id{0};
  int ver_1{0};
  int ver_2{0};
  bool operator==(const edge_& a) { return std::abs(id) == std::abs(a.id); }
  bool operator<(const edge_& a) { return std::abs(id) < std::abs(a.id); }
};

struct face_ {
  int id{0};
  std::vector<int> ver_;
  std::vector<int> edge_;
  double face_eq_d{0.0}, face_eq_a{0.0}, face_eq_b{0.0}, face_eq_c{0.0};

  bool operator==(const face_& a) { return std::abs(id) == std::abs(a.id); }
  bool operator<(const face_& a) { return std::abs(id) < std::abs(a.id); }
};

struct polyhedron_ {
  int id{0};
  std::vector<int> face_;
  bool operator==(const polyhedron_& a) { return std::abs(id) == std::abs(a.id); }
};

struct shape_ {
  std::vector<vec3r> v;
  std::vector<std::pair<size_t, size_t>> e;
  std::vector<std::vector<size_t>> f;
};

class Tess2Shape {

 public:
  Tess2Shape() {};
  virtual ~Tess2Shape() {};

  void readCommands(const char* name);

  void readTess();
  void extractShapes();
  void cutShapes();

  void writeRockableParticles(std::ostream & os);
  void writeRockableInputFile();
  void writeRockableShapeFile();

  void exec();

  // tools
  int PolyToShape(const polyhTool::poly& block, Shape& shape);
  polyhTool::poly cube(vec3r min_v, vec3r max_v);

 private:
  std::vector<vec3r> seed_positions;
  std::vector<vertex_> vertex;
  std::vector<edge_> edge;
  std::vector<face_> face;
  std::vector<polyhedron_> polyhedron;

  std::vector<Shape> Shapes;

  std::vector<std::vector<plan_>> face_plans;
  std::vector<std::vector<polyhTool::plan>> shape_plans;
  
  int number_of_cells;

  // READ IN THE COMMAND FILE
  double radius;

  std::string tessFileName{"n100-id1.tess"};
  std::string inputFileName{"input.txt"};
  std::string shapesFileName{"shapes.txt"};

  int ParticlesGroup{0};
  int ParticlesCluster{0};
  double ParticlesHomothety{1.0};

  size_t MCnstep{50000};
};

#endif /* TESS2SHAPE_HPP */
