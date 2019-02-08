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

#ifndef BREAKABLE_INTERFACE_HPP
#define BREAKABLE_INTERFACE_HPP

#include <cstdlib>     // for size_t
#include <functional>  // it holds std::less
#include <vector>

class Interaction;

class BreakableInterface {
 public:
  size_t i;  ///< ID-number of the first particle
  size_t j;  ///< ID-number of the second particle

  // Be careful not to change the order of the following parameters
  // because they are accessed in Rockable::setVariableStickParams
  // with a shift relative to 'kn'
  double kn;     ///< Embeded kn
  double kt;     ///< Embeded kt
  double kr;     ///< Embeded kr
  double fn0;    ///< Embeded fn0
  double ft0;    ///< Embeded ft0
  double mom0;   ///< Embeded mom0
  double power;  ///< Embeded power

  double dn0;    ///< distance for which fn = 0
  
  int isInner;   ///< Equals 1 if the interface is inside a cluster, 0 if it is in-between two different clusters
  std::vector<Interaction*> concernedBonds;  ///< Links to the connected points

  BreakableInterface();  // Ctor
  BreakableInterface(size_t I, size_t J);
};

namespace std {

template <>
struct less<BreakableInterface> {
  bool operator()(const BreakableInterface& lhs, const BreakableInterface& rhs) const {
    if (lhs.i < rhs.i) return true;
    if (lhs.i == rhs.i && lhs.j < rhs.j) return true;
    return false;
  }
};

template <>
struct less<BreakableInterface*> {
  bool operator()(const BreakableInterface* lhs, const BreakableInterface* rhs) const {
    if (lhs->i < rhs->i) return true;
    if (lhs->i == rhs->i && lhs->j < rhs->j) return true;
    return false;
  }
};
}  // namespace std

#endif /* end of include guard: BREAKABLE_INTERFACE_HPP */
