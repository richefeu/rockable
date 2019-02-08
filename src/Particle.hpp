// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

#include "quat.hpp"
#include "vec3.hpp"

#include "Shape.hpp"

/// A particle (which is a sphero-polyhedron)
class Particle {
 public:
  int group;    ///< A number that relates to a 'category of bodies'
  int cluster;  ///< A number that identifies the cluster to which the particle belongs

  vec3r pos;  ///< Position
  vec3r vel;  ///< Velocity
  vec3r acc;  ///< Acceleration

  quat Q;      ///< Angular position
  vec3r vrot;  ///< Angular velocity
  vec3r arot;  ///< Angular acceleration

  Shape* shape;      ///< The particle shape
  double homothety;  ///< Homothety applied to the shape
  vec3r inertia;     ///< Inertia values (same value in the diagonal)
  double mass;       ///< The particle mass

  vec3r force;   ///< Resultant force acting on the particle
  vec3r moment;  ///< Resultant moment acting on the particle

  OBB obb;  ///< precomputed OBB that fit the particle (not enlarged)

  Particle();  // Ctor

  vec3r Glob(vec3r& pos) const;
  vec3r GlobVertex(size_t v) const;
  vec3r GlobFaceVertex(size_t f, size_t v) const;
  double MinskowskiRadius() const;

  void updateObb();

  static bool VertexIsNearVertex(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax);
  static bool VertexIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax);
  static bool VertexIsNearFace(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax);
  static bool EdgeIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax);
};

#endif /* end of include guard: PARTICLE_HPP */
