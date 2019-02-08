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

#include "DataExtractor.hpp"

DataExtractor::DataExtractor() {}
DataExtractor::~DataExtractor() {}
void DataExtractor::plug(Rockable* Box) { box = Box; }
void DataExtractor::init() {}

void DataExtractor::generateHelp(std::ostream& os) {
  os << "==========================================================================================" << std::endl;
  os << "Data File: " << filename << std::endl;
  os << docString.str() << std::endl;
  for (size_t i = 0; i < columnDoc.size(); i++) {
    os << i + 1 << ":\t" << columnDoc[i] << std::endl;
  }
  os << "==========================================================================================" << std::endl
     << std::endl;
}
