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
