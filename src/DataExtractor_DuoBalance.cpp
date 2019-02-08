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
