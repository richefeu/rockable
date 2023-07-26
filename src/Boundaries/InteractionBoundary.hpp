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

#ifndef INTERACTION_BOUNDARY_HPP
#define INTERACTION_BOUNDARY_HPP

#include <cstdlib>     // for size_t
#include <functional>  // it holds std::less and std::bind

#include "Core/Particle.hpp"

//#include "Ball.hpp"

#include "vec3.hpp"

/// Data of Interaction
class InteractionBoundary {
 public:
  size_t i;     ///< Id-number of the particle
  size_t isub;  ///< id of sub-body (that can be vertex, edge or face depending on type) in Boundary i

  vec3r prev_n;  ///< Normal unit vector at the previous step
  vec3r n;       ///< Normal unit vector (from j to i)

  double dn;       ///< Normal distance (overlap = negative distance)
  double prev_dn;  ///< Backup of dn at the previous time step

  vec3r pos;  ///< Contact position
  vec3r vel;  ///< Relative velocity (j relative to i)

  double fn;  ///< Normal force (scalar value)
  vec3r ft;   ///< Tangential force (vector)
  vec3r mom;  ///< Moment at the contact point

  double damp;  ///< Precomputed vicuous damping coefficient

  virtual void deactivate();

  virtual bool update(Boundary* Pi, Particle& Pj) = 0;

  virtual InteractionBoundary& operator=(const InteractionBoundary& other);  // copy assignment operator

 protected:
  InteractionBoundary();                              // Ctor
  InteractionBoundary(const InteractionBoundary& I);  // copy Ctor
  InteractionBoundary(size_t I, size_t Isub, double Damp);
};

// This is used to sort the Interactions having the same i-particle (in
// std::set). And also for uniqueness. The lexicographic sort is organized
// by i, j, type, isub and then jsub
namespace std {
template <>
struct less<InteractionBoundary> {
  bool operator()(const InteractionBoundary& lhs, const InteractionBoundary& rhs) const {
		
#if 0
    if (lhs.i < rhs.i) return true;
    if (lhs.i == rhs.i) return (lhs.isub < rhs.isub);
    return false;
#else
    if (lhs.i != rhs.i) return lhs.i < rhs.i;
    return lhs.j < rhs.j;
#endif
		
  }
};

template <>
struct less<InteractionBoundary*> {
  bool operator()(const InteractionBoundary* lhs, const InteractionBoundary* rhs) const {

#if 0
    if (lhs->i < rhs->i) return true;
    if (lhs->i == rhs->i) return (lhs->isub < rhs->isub);
    return false;
#else
    if (lhs->i != rhs->i) return lhs->i < rhs->i;
    return lhs->j < rhs->j;
#endif
		
  }
};
}  // namespace std

#endif /* end of include guard: INTERACTION_BOUNDARY_HPP */
