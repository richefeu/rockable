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

#include "factory.hpp"
#include "kwParser.hpp"

#include "PostProcessor_ClusterGranulo.hpp"
#include "Rockable.hpp"

static Registrar<PostProcessor, ClusterGranulo> registrar("ClusterGranulo");

ClusterGranulo::ClusterGranulo() { }

void ClusterGranulo::read(std::istream& is) {
  kwParser parser;
  parser.kwMap["SievingSizes"] = __DO__(is) {
    size_t nb;
    is >> nb;
    SievingSizes.resize(nb);
    for (size_t i = 0 ; i < nb ; i++) {
      double v;
      is >> v;
      SievingSizes[i] = v;
    }
  };
  parser.parse(is);
}

void ClusterGranulo::init() {

}

void ClusterGranulo::end() {

}

void ClusterGranulo::exec() {
  char fname[256];
  sprintf(fname, "ClusterGranulo-conf%d.txt", box->iconf);
  std::ofstream file(fname);
  
  std::vector <clusterParticles> subclusters;
  box->getBrokenSubClusters(subclusters);
  if (subclusters.empty()) return;
  
  std::vector<double> NumberInClass(SievingSizes.size());
  
  for (size_t i = 0 ; i < subclusters.size() ; i++) {
    size_t sieveClass = subclusters[i].particleId.size() - 1;
    if (sieveClass >= SievingSizes.size()) sieveClass = SievingSizes.size() - 1;
    NumberInClass[sieveClass] += 1.0;
  }
  
  for (size_t i = 0 ; i < NumberInClass.size() ; i++) {
    file << SievingSizes[i] << ' ' << NumberInClass[i] / subclusters.size() << std::endl;
  }
  
}


