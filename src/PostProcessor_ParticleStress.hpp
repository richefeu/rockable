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

#ifndef POSTPROCESSOR_PARTICLESTRESS_HPP
#define POSTPROCESSOR_PARTICLESTRESS_HPP

#include "PostProcessor.hpp"

class ParticleStress : public PostProcessor {
 public:
  ParticleStress();
  void init();
  void end();
  void read(std::istream& is);

  void exec();

 private:
   double Volume;
   std::map<int, double> ConfVolumes; // associate a conf number with a volume
};

#endif /* end of include guard: POSTPROCESSOR_PARTICLESTRESS_HPP */
