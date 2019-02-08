#ifndef BODYFORCE_ATTRACTINGPOINT_HPP
#define BODYFORCE_ATTRACTINGPOINT_HPP

#include "BodyForce.hpp"

class AttractingPoint : public BodyForce {
public:
  
  AttractingPoint();
  
  void read(std::istream& is);
  void write(std::ostream& os);
  void getForceAndMoment(size_t ibody, vec3r & force, vec3r & moment);
  
private:
  vec3r point;
  double acceleration;
};

#endif /* end of include guard: BODYFORCE_ATTRACTINGPOINT_HPP */
