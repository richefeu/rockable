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

#include "DataExtractor_TrackDamage.hpp"
#include "Core/Rockable.hpp"
#include <iostream>

TrackDamage::TrackDamage() : inAreaInterfaces (0) {}

void TrackDamage::read(std::istream& is) {
  is >>  filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

  // Documentation of the output file
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("Damage");

}

void TrackDamage::init() {
  // Compute total initial area of all the interfaces
  for (auto& interfaceSet : box->Interfaces) {   
    for (auto& interface : interfaceSet) {
      inAreaInterfaces += interface.area; 
    }
  }
  std::cout << "Total initial interface area: " << inAreaInterfaces << "\n";
}



void TrackDamage::exec() {}

void TrackDamage::record() {
  double D = box->totalBrokenArea / inAreaInterfaces;
  recordFile << box->t << ' ' << D << '\n';
}

void TrackDamage::end() {}
