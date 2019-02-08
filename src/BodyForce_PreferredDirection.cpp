#include "factory.hpp"

#include "Rockable.hpp"
#include "BodyForce_PreferredDirection.hpp"

static Registrar<BodyForce, PreferredDirection> registrar("PreferredDirection");

PreferredDirection::PreferredDirection() { }

void PreferredDirection::read(std::istream& is) {
  is >> axisBody >> axis >> momentMax;
  axisBody.normalize();
  axis.normalize();
  Kr = 2.0 * momentMax / M_PI;
}

void PreferredDirection::write(std::ostream& os) {
  os << "PreferredDirection " << axisBody << ' ' << axis << ' ' << momentMax << '\n';
}

void PreferredDirection::getForceAndMoment(size_t ibody, vec3r & /*force*/, vec3r & moment) {

  vec3r axisB = box->Particles[ibody].Q * axisBody;
  axisB.normalize();
  vec3r u = cross(axis, axisB);
  u.normalize();
  double cosa = axisB * axis;
  double angle = acos(cosa);
  double sgn = -1.0;
  if (cosa < 0.0) {
    sgn = 1.0;
    angle = M_PI - angle;
  }
  
  //cout << "u  = " << u << '\n';
  //cout << "cos a = " << cosa << '\n';
  //cout << "angle = " << angle *180./M_PI << '\n';
  //cout << "sgn  = " << sgn << '\n';
  //cout << "sgn X u X angle = " << sgn*u*angle *180./M_PI << '\n';
  
  moment = Kr * sgn * u * angle;
}
