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

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <memory>

#include "quat.hpp"
#include "vec3.hpp"

#include "Shape.hpp"

struct BeemanMoreData {
  vec3r accPrevious;
  vec3r arotPrevious;
  vec3r velPrevious;
  vec3r vrotPrevious;
  vec3r accCurrent;
  vec3r arotCurrent;
};

struct RK4MoreData {
  vec3r pos0;
  quat Q0;
  vec3r vel0;
  vec3r vrot0;
  vec3r k1acc;
  vec3r k2acc;
  vec3r k3acc;
  vec3r k4acc;
  vec3r k1arot;
  vec3r k2arot;
  vec3r k3arot;
  vec3r k4arot;
};

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
  
  std::shared_ptr<BeemanMoreData> beemanData;
  std::shared_ptr<RK4MoreData> RK4Data;

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

  static bool VertexIsNearVertex(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax, const vec3r& branchPerioCorr);
  static bool VertexIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax, const vec3r& branchPerioCorr);
  static bool VertexIsNearFace(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax, const vec3r& branchPerioCorr);
  static bool EdgeIsNearEdge(Particle& Pi, Particle& Pj, size_t isub, size_t jsub, double dmax, const vec3r& branchPerioCorr);
};

#endif /* end of include guard: PARTICLE_HPP */
