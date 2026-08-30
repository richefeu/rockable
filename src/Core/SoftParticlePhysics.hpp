#pragma once

#include "mat9.hpp"

class Particle;

class SoftParticlePhysics {

  Particle* owner;
  mat9r uniformTransformation;
  mat9r stress;

  SoftParticlePhysics(Particle* t_owner);  // Ctor
};
