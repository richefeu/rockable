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

#ifndef CONTACTPARTNERSHIP_HPP
#define CONTACTPARTNERSHIP_HPP

#include <functional>
#include <map>
#include <string>
#include <utility>

class Rockable;
class Interaction;

class ContactPartnership {
 public:
  std::string name;
  std::function<void(Rockable& box)> update;        //< Compute the weights
  std::function<double(Interaction& I)> getWeight;  //< get the weight that will affect the stiffnesses

  ContactPartnership();
  ~ContactPartnership();

  void setModel(std::string& modelName);

 private:
  std::map<std::pair<size_t, size_t>, double> weightMap;
};

#endif /* end of include guard: CONTACTPARTNERSHIP_HPP */
