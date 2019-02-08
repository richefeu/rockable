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

#ifndef BLOCKRELEASE_HPP
#define BLOCKRELEASE_HPP

#include "freeFlight.hpp"
#include "vec3.hpp"
#include "quat.hpp"

struct BlockRelease {
  // The following parameters define the initial condition.
  // It should be equivalent to the initial conditions of the first free-flight
  vec3r pos;
  vec3r vel;
  quat Q;
  vec3r vrot;

  // A series of free flights
  std::vector<freeFlight> freeFlights;
};

#endif /* end of include guard: BLOCKRELEASE_HPP */
