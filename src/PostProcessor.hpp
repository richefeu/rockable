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

#ifndef POSTPROCESSOR_HPP_E4213D85
#define POSTPROCESSOR_HPP_E4213D85

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Rockable;

class PostProcessor {
 public:
  Rockable* box; ///< The Rockable instance from which data will be extracted

  virtual void plug(Rockable* Box);
  virtual void init();
  virtual void end();
  
  virtual void read(std::istream& is) = 0;
  virtual void exec() = 0;

  virtual ~PostProcessor();  // virtual Dtor

 protected:
  PostProcessor();
};

#endif /* end of include guard: POSTPROCESSOR_HPP_E4213D85 */
