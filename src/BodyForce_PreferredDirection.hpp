#ifndef BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB
#define BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB

#include "BodyForce.hpp"

class PreferredDirection : public BodyForce {
public:
  
  PreferredDirection();
  
  void read(std::istream& is);
  void write(std::ostream& os);
  void getForceAndMoment(size_t ibody, vec3r & force, vec3r & moment);
  
private:
  vec3r axisBody;
  vec3r axis;
  double momentMax;
  
  double Kr;
};

#endif /* end of include guard: BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB */
