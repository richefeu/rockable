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

#ifndef PATATRAC_HPP
#define PATATRAC_HPP

#include "Rockable.hpp"
#include "freeFlight.hpp"
#include "BlockRelease.hpp"

class Patatrac {
 public:
  Rockable box;
  std::vector<BlockRelease> releases;
  double epsilon_dn; // -> private
  size_t iblock; // -> private
  double distCol; // surrounding radius

  Patatrac();  // Ctor
  void loadDropConfigsFromFile(const char* name);
  void saveTrajectories(const char* name);
  double getDnMin();
  int getCollisionTime(freeFlight& F);
  int getCollisionTime_v2(freeFlight& F);
  void refineCollisionTime(freeFlight& F);
  int computeCollision(freeFlight& F);
  void dropAll();
  void dropAll_v2();
};

#endif /* end of include guard: PATATRAC_HPP */
