// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

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
