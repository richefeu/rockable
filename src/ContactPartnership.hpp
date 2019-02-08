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
