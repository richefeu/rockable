//  Copyright or © or Copr. Rockable
//  
//  vincent.richefeu@3sr-grenoble.fr
//  
//  This software is a computer program whose purpose is 
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces. 
//  It is developed for an ACADEMIC USAGE
//  
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use, 
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info". 
//  
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability. 
//  
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or 
//  data to be ensured and,  more generally, to use and operate it in the 
//  same conditions as regards security. 
//  
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#ifndef DRIVING_SYSTEM_HPP
#define DRIVING_SYSTEM_HPP

#include <fstream>
#include <functional>
#include <utility>
#include <vector>

#include "vec3.hpp"
#include "mat9.hpp"

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
#define _xrot_Mom_ 9
#define _yrot_Mom_ 10
#define _zrot_Mom_ 11

// values larger or equal to 100 are imposed vector (rather than single value)
#define _xyzrot_Vel_ 100
#define _xyzrot_Mom_ 101

/// @brief Control of ONE degree of freedom
struct Control {
  int type;      // type of control (_x_Vel_, etc.)
  size_t i;      // identifier number of the controlled particle
  double value;  // imposed value (translation/rotation velocity or force/moment)
  vec3r vec_value; // in case we impose a velocity or force vector
};

struct PeriodicCellControl {
  mat9b Drive; ///< Driving mode. Can be ForceDriven or VelocityDriven
  mat9r Sig;   ///< Imposed external stress
  mat9r v;     ///< Imposed velocities
};

/// @brief This class hold the controls (force or velocity on the 'nContr' first
/// bodies)
class DrivingSystem {
 public:
  std::function<void(Rockable& box)> ServoFunction;  ///< This lamda function is
                                                     ///< used for 'smart' driving conditions

  std::vector<Control> controls;
  PeriodicCellControl cellControl;

  DrivingSystem();  // Ctor

  void read(bool allow_warn = true);
};

#endif /* end of include guard: DRIVING_SYSTEM_HPP */
