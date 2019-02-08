//        Rockable, 3D-DEM with sphero-polyhedra
//        Copyright (C) 2016-2019  <vincent.richefeu@3sr-grenoble.fr>
//        
//        This program is free software: you can redistribute it and/or modify
//        it under the terms of the GNU General Public License as published by
//        the Free Software Foundation, either version 3 of the License, or
//        (at your option) any later version.
//        
//        This program is distributed in the hope that it will be useful,
//        but WITHOUT ANY WARRANTY; without even the implied warranty of
//        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//        GNU General Public License for more details.
//        
//        You should have received a copy of the GNU General Public License
//        along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
