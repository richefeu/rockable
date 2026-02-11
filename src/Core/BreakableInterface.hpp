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

#ifndef BREAKABLE_INTERFACE_HPP
#define BREAKABLE_INTERFACE_HPP

#include <cstdlib>     // for size_t
#include <functional>  // it holds std::less
#include <vector>

class Interaction;

const int breakModel_YieldSurface = 0;
const int breakModel_Gc = 1;

class BreakableInterface {
 public:
  size_t i{0};  ///< ID-number of the first particle
  size_t j{0};  ///< ID-number of the second particle

  // Be careful not to change the order of the following parameters
  // because they are accessed in Rockable::setVariableStickParams
  // with a shift relative to 'kn'
  double kn{0.0};     ///< Embeded kn
  double kt{0.0};     ///< Embeded kt
  double kr{0.0};     ///< Embeded kr
  double fn0{0.0};    ///< Embeded fn0
  double ft0{0.0};    ///< Embeded ft0
  double mom0{0.0};   ///< Embeded mom0
  double power{0.0};  ///< Embeded power
  double Gc{0.0};     ///< Embeded surface energy

  int breakModel{breakModel_Gc};

  double dn0{0.0};   ///< distance for which fn = 0
  double area{0.0};  ///< ___
  double En{0.0};    ///< ___
  double Et{0.0};    ///< ___

  int isInner{-1};  ///< Equals 1 if the interface is inside a cluster, 0 if it is in-between two different clusters
  std::vector<Interaction*> concernedBonds;  ///< Links to the connected points

  BreakableInterface();  // Ctor
  BreakableInterface(size_t I, size_t J);
};

namespace std {

template <>
struct less<BreakableInterface> {
  bool operator()(const BreakableInterface& lhs, const BreakableInterface& rhs) const {
    if (lhs.i < rhs.i) {
      return true;
    }
    if (lhs.i == rhs.i && lhs.j < rhs.j) {
      return true;
    }
    return false;
  }
};

template <>
struct less<BreakableInterface*> {
  bool operator()(const BreakableInterface* lhs, const BreakableInterface* rhs) const {
    if (lhs->i < rhs->i) {
      return true;
    }
    if (lhs->i == rhs->i && lhs->j < rhs->j) {
      return true;
    }
    return false;
  }
};
}  // namespace std

#endif /* end of include guard: BREAKABLE_INTERFACE_HPP */
