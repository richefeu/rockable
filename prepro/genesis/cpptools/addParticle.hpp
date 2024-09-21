#ifndef ADD_PARTICLE_HPP
#define ADD_PARTICLE_HPP

// RECALL:
// globalTransformation and individualParticleRotation are global variables
// this .h file is included in the generator.cpp 

#include <cmath>
#include <iostream>

#include "quat.hpp"
#include "vec3.hpp"
#include "transformation.hpp"

extern Transformation<double> globalTransformation;
extern quat individualParticleRotation;

void addParticle(std::ostream& os, const char* name, int group, int cluster, double homothety, vec3r& position,
                 quat& angularPosition) {
  using namespace std;

  os << name << ' ' << group << ' ' << cluster << ' ' << homothety << "  " << position << "  0 0 0  0 0 0   "
     << angularPosition << "  0 0 0  0 0 0\n";
}

void addParticle(std::ostream& os, const char* name, int group, int cluster, double homothety, vec3r& position) {
  using namespace std;

  os << name << ' ' << group << ' ' << cluster << ' ' << homothety << "  " << position
     << "  0 0 0  0 0 0   1 0 0 0  0 0 0  0 0 0\n";
}

#endif /* end of include guard: ADD_PARTICLE_HPP */
