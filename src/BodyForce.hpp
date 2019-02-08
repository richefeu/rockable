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
#ifndef BODYFORCE_HPP_E253559A
#define BODYFORCE_HPP_E253559A

#include <fstream>
#include "vec3.hpp"

class Rockable;

// 
class BodyForce {
public:
  Rockable* box; 
  
  virtual ~BodyForce();
  virtual void plug(Rockable* Box);
  virtual void read(std::istream& is) = 0;
  virtual void write(std::ostream& os) = 0;
  virtual void getForceAndMoment(size_t ibody, vec3r & force, vec3r & moment) = 0;
  
protected:
  BodyForce();
};

#endif /* end of include guard: BODYFORCE_HPP_E253559A */