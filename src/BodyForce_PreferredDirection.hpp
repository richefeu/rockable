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

#ifndef BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB
#define BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB

#include "BodyForce.hpp"

class PreferredDirection : public BodyForce {
public:
  
  PreferredDirection();
  
  void read(std::istream& is);
  void write(std::ostream& os);
  void getForceAndMoment(size_t ibody, vec3r & force, vec3r & moment);
  
private:
  vec3r axisBody;
  vec3r axis;
  double momentMax;
  
  double Kr;
};

#endif /* end of include guard: BODYFORCE_PREFERREDDIRECTION_HPP_021F7DFB */
