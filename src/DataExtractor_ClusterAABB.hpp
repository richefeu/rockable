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

#ifndef CLUSTERAABB_HPP_7FF11E00
#define CLUSTERAABB_HPP_7FF11E00

#include <fstream>

#include "DataExtractor.hpp"

class ClusterAABB : public DataExtractor {
 public:
  ClusterAABB();
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  int icluster;
};

#endif /* end of include guard: CLUSTERAABB_HPP_7FF11E00 */
