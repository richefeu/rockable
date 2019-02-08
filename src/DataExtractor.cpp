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

#include "DataExtractor.hpp"

DataExtractor::DataExtractor() {}
DataExtractor::~DataExtractor() {}
void DataExtractor::plug(Rockable* Box) { box = Box; }
void DataExtractor::init() {}

void DataExtractor::generateHelp(std::ostream& os) {
  os << "==========================================================================================" << std::endl;
  os << "Data File: " << filename << std::endl;
  os << docString.str() << std::endl;
  for (size_t i = 0; i < columnDoc.size(); i++) {
    os << i + 1 << ":\t" << columnDoc[i] << std::endl;
  }
  os << "==========================================================================================" << std::endl
     << std::endl;
}
