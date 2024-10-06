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

#ifndef ADD_PARTICLE_HPP
#define ADD_PARTICLE_HPP

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp

#include <cmath>
#include <iostream>

#include "quat.hpp"
#include "transformation.hpp"
#include "vec3.hpp"

void addParticle(std::ostream& os, const char* name, int group, int cluster, double homothety, vec3r& position,
                 quat& angularPosition) {
  using namespace std;
  vec3r pos = position;
  globalTransformation.apply(pos);
  os << name << ' ' << group << ' ' << cluster << ' ' << homothety << "  " << pos << "  0 0 0  0 0 0   "
     << angularPosition * individualParticleRotation << "  0 0 0  0 0 0\n";
}

// DEPRECATED (TO REMOVE...)
/*
void addParticle(std::ostream& os, const char* name, int group, int cluster, double homothety, vec3r& position) {
  using namespace std;

  vec3r pos = position;
  globalTransformation.apply(pos);
  os << name << ' ' << group << ' ' << cluster << ' ' << homothety << "  " << pos
     << "  0 0 0  0 0 0   1 0 0 0  0 0 0  0 0 0\n";
}
*/

#endif /* end of include guard: ADD_PARTICLE_HPP */
