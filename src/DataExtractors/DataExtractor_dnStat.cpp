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

#include "DataExtractor_dnStat.hpp"
#include "Core/Rockable.hpp"

dnStat::dnStat() {}

void dnStat::read(std::istream& is) {
  is >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

  // Documentation of the output file
  docString << "nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("dnMin");
  columnDoc.push_back("dnMax");
  columnDoc.push_back("dnMean");
  columnDoc.push_back("dnMeanNeg");
  columnDoc.push_back("dnMeanPos");
  columnDoc.push_back("nbNeg");
  columnDoc.push_back("nbPos");
  columnDoc.push_back("pos.x of dnMin");
  columnDoc.push_back("pos.y of dnMin");
  columnDoc.push_back("pos.z of dnMin");
}

void dnStat::init() {}

void dnStat::exec() {}

void dnStat::record() {
  double dnMin = 1e20;
  double dnMax = -1e20;
  double dnMeanNeg = 0.0;
  double dnMeanPos = 0.0;
  double nbNeg = 0.0;
  double nbPos = 0.0;
  double dnMean = 0.0;
  vec3r posMin;
  for (size_t i = 0; i < box->activeInteractions.size(); i++) {
    double dn = box->activeInteractions[i]->dn;
    dnMean += dn;
    if (dn < 0.0) {
      dnMeanNeg += dn;
      nbNeg += 1.0;
    } else {
      dnMeanPos += dn;
      nbPos += 1.0;
    }
    if (dn < dnMin) {
      dnMin = dn;
      posMin = box->activeInteractions[i]->pos;
    }
    if (dn > dnMax) {
      dnMax = dn;
    }
  }
  if (nbNeg > 0.0) dnMeanNeg /= nbNeg;
  if (nbPos > 0.0) dnMeanPos /= nbPos;
  if (nbNeg > 0.0 || nbPos > 0.0) dnMean /= (nbNeg + nbPos);

  recordFile << box->t << ' ' << dnMin << ' ' << dnMax << ' ' << dnMean << ' ' << dnMeanNeg << ' ' << dnMeanPos << ' '
             << (int)nbNeg << ' ' << (int)nbPos << ' ' << posMin << std::endl
             << std::flush;
}

void dnStat::end() {}
