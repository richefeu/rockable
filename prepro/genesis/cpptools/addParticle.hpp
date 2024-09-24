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
