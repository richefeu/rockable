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

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <fstream>
#include <functional>

#include "OBB.hpp"
#include "vec3.hpp"

// A Boundary is a kind of particle with special shape.
// It is concave or flat so that only interactions with node spheres are considered. 
// This is why its structure is so similar to that of a particle.
class Boundary {
 public:
  int group;            ///< A number that relates to a 'category of bodies'
  double mass_masstot;  ///< ???
  double masstot;       ///< The particle mass
  vec3r pos;            ///< Position
  double force;         ///< Resultant force acting on the particle

  double vel;
  double acc;
  double sig;

  vec3r vrot;
  vec3r arot;

  double damp;

  unsigned int nbP;

  virtual ~Boundary();
  virtual void read(std::istream& is) = 0;
  virtual void write(std::ostream& os) = 0;

  virtual bool SphereIsNear(vec3r& position, double radius, double dmax) = 0;
  virtual bool ObbIsNear(OBB& Pj, double dmax) = 0;

 protected:
  Boundary();
};

#endif /* end of include guard: BOUNDARY_HPP */