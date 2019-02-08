// Copyright (C) Patatrac <vincent.richefeu@3sr-grenoble.fr>
// 
// This file is part of mbox.
// 
// Patatrac can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note 
// Without a license, the code is copyrighted by default. 
// People can read the code, but they have no legal right to use it. 
// To use the code, you must contact the author directly and ask permission.

#ifndef PATATRAC_HPP
#define PATATRAC_HPP

#include "Rockable.hpp"
#include "freeFlight.hpp"
#include "BlockRelease.hpp"

class Patatrac {
 public:
  Rockable box;
  std::vector<BlockRelease> releases;
  double epsilon_dn; // -> private
  size_t iblock; // -> private
  double distCol; // surrounding radius

  Patatrac();  // Ctor
  void loadDropConfigsFromFile(const char* name);
  void saveTrajectories(const char* name);
  double getDnMin();
  int getCollisionTime(freeFlight& F);
  int getCollisionTime_v2(freeFlight& F);
  void refineCollisionTime(freeFlight& F);
  int computeCollision(freeFlight& F);
  void dropAll();
  void dropAll_v2();
};

#endif /* end of include guard: PATATRAC_HPP */
