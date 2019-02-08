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
