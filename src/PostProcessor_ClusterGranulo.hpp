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

#ifndef POSTPROCESSOR_CLUSTERGRANULO_CPP
#define POSTPROCESSOR_CLUSTERGRANULO_CPP

#include "PostProcessor.hpp"

class ClusterGranulo : public PostProcessor {
 public:
  ClusterGranulo();
  void init();
  void end();
  void read(std::istream& is);

  void exec();

 private:
   std::vector<double> SievingSizes;
};

#endif /* end of include guard: POSTPROCESSOR_CLUSTERGRANULO_CPP */
