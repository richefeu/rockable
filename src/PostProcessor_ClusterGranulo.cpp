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


