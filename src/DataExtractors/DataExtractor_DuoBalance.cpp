//  Copyright or © or Copr. Rockable
//  
//  vincent.richefeu@3sr-grenoble.fr
//  
//  This software is a computer program whose purpose is 
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces. 
//  It is developed for an ACADEMIC USAGE
//  
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use, 
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info". 
//  
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability. 
//  
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or 
//  data to be ensured and,  more generally, to use and operate it in the 
//  same conditions as regards security. 
//  
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include <limits>

#include "factory.hpp"

#include "DataExtractor_DuoBalance.hpp"
#include "Core/Rockable.hpp"

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
