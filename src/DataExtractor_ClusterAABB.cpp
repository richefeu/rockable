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
