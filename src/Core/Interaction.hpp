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

#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include <cstdlib>     // for size_t
#include <functional>  // it holds std::less and std::bind

#include "Particle.hpp"
#include "vec3.hpp"

//====================   isub is:     jsub is:  ===========
const int vvType = 0;  // vertex       vertex
const int veType = 1;  // vertex       edge
const int vfType = 2;  // vertex       face
const int eeType = 3;  // edge         edge
//=========================================================
const std::string InteractionTypeName[] = {"vertex-vertex", "vertex-edge", "vertex-face", "edge-edge"};
//=========================================================

class BreakableInterface;

/// Data of Interaction
class Interaction {
 public:
  size_t i;     ///< ID-number of the first particle
  size_t j;     ///< ID-number of the second particle
                // !! BE CAREFUL: there's no explicit assertion that i < j
  int type;     ///< The type of contact
  size_t isub;  ///< id of sub-body (that can be vertex, edge or face depending
                ///< on type) in sphero-polyhedron i
  size_t jsub;  ///< id of sub-body (that can be vertex, edge or face depending
                ///< on type) in sphero-polyhedron j

  vec3r prev_n;  ///< Normal unit vector at the previous step
  vec3r n;       ///< Normal unit vector (from j to i)

  double dn{0.0};       ///< Normal distance (overlap = negative distance)
  double prev_dn{0.0};  ///< Backup of dn at the previous time step

  vec3r ds;

  vec3r pos;             ///< Contact position
  vec3r vel;             ///< Relative velocity (j relative to i)
  vec3r jPeriodicShift;  ///< In case of periodic cell, this the shifting of particle j

  double fn{0.0};  ///< Normal force (scalar value)
  vec3r ft;        ///< Tangential force (vector)
  vec3r mom;       ///< Moment at the contact point

  double damp{0.0};           ///< Precomputed vicuous damping coefficient
  BreakableInterface* stick{nullptr};  ///< The breakable interface (if there's no interface, then stick = nullptr)

  Interaction();                      // Ctor
  Interaction(const Interaction& I);  // copy Ctor
  Interaction(size_t I, size_t J, int Type, size_t Isub, size_t Jsub, double Damp, BreakableInterface* Stick = nullptr);
  Interaction& operator=(const Interaction& other);  // copy assignment operator
  void deactivate();

  // UpdateDispatcher is the update-dispatcher as a function of the four types
  // (vv, ve, vf, ee). Note: Krishna Kumar advised me to use lambda instead of
  // binded functions in UpdateDispatcher. This change provided a gain (redution
  // of CPU time) of about 12% for the tested example.
  static std::function<bool(Interaction&, Particle&, Particle&)> UpdateDispatcher[4];
  static std::function<bool(Interaction&, Particle&, Particle&)> UpdateDispatcherPeriodic[4];
};

// This is used to sort the Interactions having the same i-particle (in
// std::set). And also for uniqueness. The lexicographic sort is organized
// by i, j, type, isub and then jsub
namespace std {
template <>
struct less<Interaction> {
  bool operator()(const Interaction& lhs, const Interaction& rhs) const {
    if (lhs.i != rhs.i) {
      return lhs.i < rhs.i;
    }
    if (lhs.j != rhs.j) {
      return lhs.j < rhs.j;
    }
    if (lhs.type != rhs.type) {
      return lhs.type < rhs.type;
    }
    if (lhs.isub != rhs.isub) {
      return lhs.isub < rhs.isub;
    }
    return lhs.jsub < rhs.jsub;
  }
};

template <>
struct less<Interaction*> {
  bool operator()(const Interaction* lhs, const Interaction* rhs) const {
    if (lhs->i != rhs->i) {
      return lhs->i < rhs->i;
    }
    if (lhs->j != rhs->j) {
      return lhs->j < rhs->j;
    }
    if (lhs->type != rhs->type) {
      return lhs->type < rhs->type;
    }
    if (lhs->isub != rhs->isub) {
      return lhs->isub < rhs->isub;
    }
    return lhs->jsub < rhs->jsub;
  }
};
}  // namespace std

#endif /* end of include guard: INTERACTION_HPP */
