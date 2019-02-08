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

#include "AABB.hpp"
#include "factory.hpp"

#include "DataExtractor_ClusterAABB.hpp"
#include "Rockable.hpp"

static Registrar<DataExtractor, ClusterAABB> registrar("ClusterAABB");

ClusterAABB::ClusterAABB() : icluster(0) {}

void ClusterAABB::read(std::istream& is) {
  is >> icluster >> filename >> nrec;
  if (box->isInteractive() == false) recordFile.open(filename.c_str());
  nstep = std::numeric_limits<int>::max();

  // Documentation of the output file
  docString << "icluster = " << icluster << ", nrec = " << nrec;
  columnDoc.clear();
  columnDoc.push_back("Time");
  columnDoc.push_back("aabb.min.x");
  columnDoc.push_back("aabb.min.y");
  columnDoc.push_back("aabb.min.z");
  columnDoc.push_back("aabb.max.x");
  columnDoc.push_back("aabb.max.y");
  columnDoc.push_back("aabb.max.z");
}

void ClusterAABB::exec() {}

void ClusterAABB::record() {
  AABB aabb;
  size_t i0 = 0;
  for (size_t i = 0; i < box->Particles.size(); i++) {
    if (box->Particles[i].cluster == icluster) {
      aabb.set_single(box->Particles[i].GlobVertex(0));
      double radius = box->Particles[i].MinskowskiRadius();
      aabb.enlarge(radius);
      for (size_t v = 0; v < box->Particles[i].shape->vertex.size(); v++) {
        AABB boxi;
        boxi.set_single(box->Particles[i].GlobVertex(v));
        boxi.enlarge(radius);
        aabb.enlarge(boxi);
      }
      i0 = i;
      break;
    }
  }

  for (size_t i = i0; i < box->Particles.size(); i++) {
    if (box->Particles[i].cluster == icluster) {
      double radius = box->Particles[i].MinskowskiRadius();
      for (size_t v = 0; v < box->Particles[i].shape->vertex.size(); v++) {

        AABB boxi;
        boxi.set_single(box->Particles[i].GlobVertex(v));
        boxi.enlarge(radius);

        aabb.enlarge(boxi);
      }
    }
  }

  recordFile << box->t << ' ' << aabb.min << ' ' << aabb.max << std::endl << std::flush;
}

void ClusterAABB::end() {}
