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

#include "factory.hpp"
#include "kwParser.hpp"

#include "PostProcessor_ClusterGranulo.hpp"
#include "Rockable.hpp"
#include "processingTool_getBrokenSubClusters.hpp"

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
  getBrokenSubClusters(box,subclusters);
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


