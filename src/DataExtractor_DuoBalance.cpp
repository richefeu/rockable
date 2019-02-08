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

#include <limits>

#include "factory.hpp"

#include "DataExtractor_DuoBalance.hpp"
#include "Rockable.hpp"

static Registrar<DataExtractor, DuoBalance> registrar("DuoBalance");

DuoBalance::DuoBalance() {}

void DuoBalance::read(std::istream& is) {
  is >> i >> j >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

  // Documentation of the output file
  docString << "first particle = " << i << ", second particle = " << j << ", nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("nvv");
  columnDoc.push_back("nve");
  columnDoc.push_back("nvf");
  columnDoc.push_back("nee");
  columnDoc.push_back("(nvv + nve + nvf + nee) x weight");
}

void DuoBalance::init() {}

void DuoBalance::exec() {}

void DuoBalance::record() {

  size_t nvv = 0;
  size_t nve = 0;
  size_t nvf = 0;
  size_t nee = 0;

  /*
  for (auto it = box->Interactions[i].begin(); it != box->Interactions[i].end(); ++it) {
    Interaction* I = const_cast<Interaction*>(std::addressof(*it));
    if (( ((I->i == i) && (I->j == j)) || ((I->i == j) && (I->j == i)) ) 
        && I->dn < 0.0) {
      if (I->type == vvType)
        nvv++;
      else if (I->type == veType)
        nve++;
      else if (I->type == vfType)
        nvf++;
      else if (I->type == eeType)
        nee++;
    }
  }
  */
  
  for (auto it = box->Interactions[i].begin(); it != box->Interactions[i].end(); ++it) {
    Interaction* I = const_cast<Interaction*>(std::addressof(*it));
    if ((I->i == j || I->j == j) && I->dn < 0.0) {
      if (I->type == vvType)
        nvv++;
      else if (I->type == veType)
        nve++;
      else if (I->type == vfType)
        nvf++;
      else if (I->type == eeType)
        nee++;
    }
  }

  for (auto it = box->Interactions[j].begin(); it != box->Interactions[j].end(); ++it) {
    Interaction* I = const_cast<Interaction*>(std::addressof(*it));
    if ((I->i == i || I->j == i) && I->dn < 0.0) {
      if (I->type == vvType)
        nvv++;
      else if (I->type == veType)
        nve++;
      else if (I->type == vfType)
        nvf++;
      else if (I->type == eeType)
        nee++;
    }
  }
  
  recordFile << box->t << "   " << nvv << ' ' << nve << ' ' << nvf << ' ' << nee << "\n";
  
  /*
  if (box->ctcPartnership.getWeight != nullptr) {
    for (auto it = box->Interactions[i].begin(); it != box->Interactions[i].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if (( ((I->i == i) && (I->j == j)) || ((I->i == j) && (I->j == i)) ) 
          && I->dn < 0.0) {
        recordFile << InteractionTypeName[I->type] << ' ' 
          << I->i << ' ' << I->j << ' ' << I->isub << ' ' << I->jsub << ' '
          << box->ctcPartnership.getWeight(*I) << ' ' 
          << I->dn << ' ' << I->fn << ' ' << I->ft << '\n';
      }
    }
  }
  */
  
  if (box->ctcPartnership.getWeight != nullptr) {
    for (auto it = box->Interactions[i].begin(); it != box->Interactions[i].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if ((I->i == j || I->j == j) && I->dn < 0.0) {
        recordFile << InteractionTypeName[I->type] << ' ' 
          << I->i << ' ' << I->j << ' ' << I->isub << ' ' << I->jsub << ' '
          << box->ctcPartnership.getWeight(*I) << ' ' 
          << I->dn << ' ' << I->fn << ' ' << I->ft << '\n';
      }
    }

    for (auto it = box->Interactions[j].begin(); it != box->Interactions[j].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));
      if ((I->i == i || I->j == i) && I->dn < 0.0) {
        recordFile << InteractionTypeName[I->type] << ' ' 
          << I->i << ' ' << I->j << ' ' << I->isub << ' ' << I->jsub << ' '
          << box->ctcPartnership.getWeight(*I) << ' ' 
          << I->dn << ' ' << I->fn << ' ' << I->ft << '\n';
      }
    }
  }
  
  recordFile << '\n' << std::flush;
}

void DuoBalance::end() {}
