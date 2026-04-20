#pragma once

#include <cmath>
#include <algorithm>

#include "vec3.hpp"

struct Slicing {
  vec3<float> reducedPos;   // point sur le plan mapper sur la boite de la scene (coordonnées entre 0 et 1)
  float azimuthDeg{0.0};    // Azimuth angle (horizontal plane, in degrees, 0-360)
  float elevationDeg{0.0};  // Elevation angle (vertical tilt, in degrees, -90, 90)
  float reducedHalfWidth{0.1};

  vec3r P;               // point sur le plan
  vec3r N;               // normale unitaire au plan
  float halfWidth{0.0};  // demi-épaisseur de la slice

  Slicing();
  void updatePos(vec3r& min, vec3r& max);
  void updateN();
  void updateWidth(vec3r& min, vec3r& max);
  bool inSlice(vec3r& pos);
};
