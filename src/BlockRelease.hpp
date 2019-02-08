#ifndef BLOCKRELEASE_HPP
#define BLOCKRELEASE_HPP

#include "freeFlight.hpp"
#include "vec3.hpp"
#include "quat.hpp"

struct BlockRelease {
  // The following parameters define the initial condition.
  // It should be equivalent to the initial conditions of the first free-flight
  vec3r pos;
  vec3r vel;
  quat Q;
  vec3r vrot;

  // A series of free flights
  std::vector<freeFlight> freeFlights;
};

#endif /* end of include guard: BLOCKRELEASE_HPP */
