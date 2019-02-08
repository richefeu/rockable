//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
