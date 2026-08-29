#pragma once

#include "AABB.hpp"
#include "vec3.hpp"

class Rockable;

struct Probing {
  AABB aabb;
  vec3r reducedPos;
  vec3r reducedSizes;
  
  Probing();
  void extractData(Rockable & box, const std::string & which_data);
};
