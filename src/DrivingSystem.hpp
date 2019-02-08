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

#ifndef DRIVING_SYSTEM_HPP
#define DRIVING_SYSTEM_HPP

#include <fstream>
#include <functional>
#include <utility>
#include <vector>

class Rockable;

const bool ForceDriven = true;
const bool VelocityDriven = false;

#define _x_Vel_ 0
#define _y_Vel_ 1
#define _z_Vel_ 2
#define _xrot_Vel_ 3
#define _yrot_Vel_ 4
#define _zrot_Vel_ 5
#define _x_For_ 6
#define _y_For_ 7
#define _z_For_ 8

// Not implemented yet (actually, not yet necessary)
#define _xrot_Mom_ 9
#define _yrot_Mom_ 10
#define _zrot_Mom_ 11

/// @brief Control of ONE degree of freedom
struct Control {
  int type;      // type of control (_x_Vel, etc.)
  size_t i;      // identifier number of the controlled particle
  double value;  // imposed value (translation/rotation velocity or force/moment)
};

/// @brief This class hold the controls (force or velocity on the 'nContr' first
/// bodies)
class DrivingSystem {
 public:
  std::function<void(Rockable& box)> ServoFunction;  ///< This lamda function is
                                                     ///< used for 'smart' driving conditions

  std::vector<Control> controls;

  DrivingSystem();  // Ctor

  void read(bool warn = true);
};

#endif /* end of include guard: DRIVING_SYSTEM_HPP */
