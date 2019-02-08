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

#ifndef DATAEXTRACTOR_HPP_E4213D85
#define DATAEXTRACTOR_HPP_E4213D85

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Rockable;

class DataExtractor {
 public:
  Rockable* box;             ///< The Rockable instance from which data will be extracted
  std::string filename;      ///< Name of a file to store extracted data
  std::ofstream recordFile;  ///< stream to a file for storing extracted data

  int nstep;  ///< Call-period for exec
  int nrec;   ///< Call-period for record

  virtual void plug(Rockable* Box);
  virtual void init();

  virtual void read(std::istream& is) = 0;

  virtual void exec() = 0;    // In case there are some computations to do with nstep periodicity
  virtual void record() = 0;  // Record with nrec periodicity
  virtual void end() = 0;     // Called after all steps have been done

  virtual ~DataExtractor();  // virtual Dtor

  virtual void generateHelp(std::ostream& os);

  std::stringstream docString;
  std::vector<std::string> columnDoc;

 protected:
  DataExtractor();
};

#endif /* end of include guard: DATAEXTRACTOR_HPP_E4213D85 */
