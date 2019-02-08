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

#ifndef POSTPROCESSOR_HPP_E4213D85
#define POSTPROCESSOR_HPP_E4213D85

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Rockable;

class PostProcessor {
 public:
  Rockable* box; ///< The Rockable instance from which data will be extracted

  virtual void plug(Rockable* Box);
  virtual void init();
  virtual void end();
  
  virtual void read(std::istream& is) = 0;
  virtual void exec() = 0;

  virtual ~PostProcessor();  // virtual Dtor

 protected:
  PostProcessor();
};

#endif /* end of include guard: POSTPROCESSOR_HPP_E4213D85 */
