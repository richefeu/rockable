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

#ifndef CLUSTERAABB_HPP_7FF11E00
#define CLUSTERAABB_HPP_7FF11E00

#include <fstream>

#include "DataExtractor.hpp"

class ClusterAABB : public DataExtractor {
 public:
  ClusterAABB();
  void read(std::istream& is);

  void exec();
  void record();
  void end();

 private:
  int icluster;
};

#endif /* end of include guard: CLUSTERAABB_HPP_7FF11E00 */
