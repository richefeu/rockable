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

#ifndef DATAEXTRACTOR_DUOBALANCE_HPP
#define DATAEXTRACTOR_DUOBALANCE_HPP

#include "DataExtractor.hpp"

class DuoBalance : public DataExtractor {
 public:
  DuoBalance();
  void init();
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  size_t i, j;
};

#endif /* end of include guard: DATAEXTRACTOR_DUOBALANCE_HPP */
