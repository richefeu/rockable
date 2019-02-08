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

#ifndef DATAEXTRACTOR_HPP_E4213D85
#define DATAEXTRACTOR_HPP_E4213D85

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Rockable;

class DataExtractor {
 public:
  Rockable* box;             ///< The Rockable instance from which data will be extracted
  std::string filename;      ///< Name of a file to store extracted data
  std::ofstream recordFile;  ///< stream to a file for storing extracted data

  int nstep;  ///< Call-period for exec
  int nrec;   ///< Call-period for record

  virtual void plug(Rockable* Box);
  virtual void init();

  virtual void read(std::istream& is) = 0;

  virtual void exec() = 0;    // In case there are some computations to do with nstep periodicity
  virtual void record() = 0;  // Record with nrec periodicity
  virtual void end() = 0;     // Called after all steps have been done

  virtual ~DataExtractor();  // virtual Dtor

  virtual void generateHelp(std::ostream& os);

  std::stringstream docString;
  std::vector<std::string> columnDoc;

 protected:
  DataExtractor();
};

#endif /* end of include guard: DATAEXTRACTOR_HPP_E4213D85 */
