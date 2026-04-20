#include "Slicing.hpp"

Slicing::Slicing() {
  reducedPos.set(0.5, 0.5, 0.5);
  // N.set(0.0, 0.0, 1.0);
  updateN();
  std::cout << "N = " << N << std::endl;
}

void Slicing::updatePos(vec3r& min, vec3r& max) {
  vec3r D = max - min;
  P.x = min.x + reducedPos.x * D.x;
  P.y = min.y + reducedPos.y * D.y;
  P.z = min.z + reducedPos.z * D.z;
}

void Slicing::updateN() {
  // Conversion en radians
  float azimuthRad = azimuthDeg * M_PI / 180.0f;
  float elevationRad = elevationDeg * M_PI / 180.0f;

  // Vecteur initial : (1, 0, 0)
  N.x = cos(elevationRad) * sin(azimuthRad);
  N.y = cos(elevationRad) * cos(azimuthRad);
  N.z = sin(elevationRad);
}

void Slicing::updateWidth(vec3r& min, vec3r& max) {
  vec3r D = max - min;
  double L = 0.5 * std::min({D.x, D.y, D.z});
  halfWidth = reducedHalfWidth * L;
}

bool Slicing::inSlice(vec3r& pos) { return (fabs((pos - P) * N) <= halfWidth); }