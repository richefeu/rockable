#include "factory.hpp"

#include "Rockable.hpp"
#include "BodyForce_AttractingPoint.hpp"

static Registrar<BodyForce, AttractingPoint> registrar("AttractingPoint");

AttractingPoint::AttractingPoint() { }

void AttractingPoint::read(std::istream& is) {
  is >> point >> acceleration;
}

void AttractingPoint::write(std::ostream& os) {
  os << "AttractingPoint " << point << ' ' << acceleration << '\n';
}

void AttractingPoint::getForceAndMoment(size_t ibody, vec3r & force, vec3r & /*moment*/) {
  vec3r direction = point - box->Particles[ibody].pos;
  double l2 = norm2(direction);
  
  if (l2 > 0.0) {
    direction *= 1.0 / sqrt(l2); // normalize
    force = box->Particles[ibody].mass * acceleration * direction;
  } 
}
