// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.
#ifndef BODYFORCE_HPP_E253559A
#define BODYFORCE_HPP_E253559A

#include <fstream>
#include "vec3.hpp"

class Rockable;

// 
class BodyForce {
public:
  Rockable* box; 
  
  virtual ~BodyForce();
  virtual void plug(Rockable* Box);
  virtual void read(std::istream& is) = 0;
  virtual void write(std::ostream& os) = 0;
  virtual void getForceAndMoment(size_t ibody, vec3r & force, vec3r & moment) = 0;
  
protected:
  BodyForce();
};

#endif /* end of include guard: BODYFORCE_HPP_E253559A */