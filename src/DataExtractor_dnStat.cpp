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

#include "DataExtractor_dnStat.hpp"
#include "Rockable.hpp"

static Registrar<DataExtractor, dnStat> registrar("dnStat");

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
