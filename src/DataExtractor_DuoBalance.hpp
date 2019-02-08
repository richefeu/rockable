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

#ifndef DATAEXTRACTOR_DUOBALANCE_HPP
#define DATAEXTRACTOR_DUOBALANCE_HPP

#include "DataExtractor.hpp"

class DuoBalance : public DataExtractor {
 public:
  DuoBalance();
  void init();
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  size_t i, j;
};

#endif /* end of include guard: DATAEXTRACTOR_DUOBALANCE_HPP */
