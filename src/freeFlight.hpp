#ifndef FREEFLIGHT_HPP
#define FREEFLIGHT_HPP

#include "vec3.hpp"
#include "quat.hpp"

// A free flight is an initial configuration (pos, vel, Q, vrot),
// a given time ti and a duration.
struct freeFlight {
  double ti;  // Initial time of the free-flight
  double duration;

  // kinematic-values at initial time ti
  vec3r pos;
  vec3r vel;
  quat Q;
  vec3r vrot;
};

#endif /* end of include guard: FREEFLIGHT_HPP */
