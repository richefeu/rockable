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

#include "Core/Rockable.hpp"
#include "PreproCommand_copyParamsToInterfaces.hpp"

copyParamsToInterfaces::copyParamsToInterfaces() { }

void copyParamsToInterfaces::addCommand() {  
  box->parser.kwMap["copyParamsToInterfaces"] = [this](std::istream& conf) {
    conf >> this->isInnerStr;
    exec();
  };
}

/**
   @attention  We suppose that the Damp parameter has already been set
               in the corresponding Interactions
*/
void copyParamsToInterfaces::exec() {
  int isInner = 0;
  if (isInnerStr == "inner") isInner = 1;

  for (size_t i = 0; i < box->Interfaces.size(); i++) {
    for (auto it = box->Interfaces[i].begin(); it != box->Interfaces[i].end(); ++it) {
      int icluster = box->Particles[it->i].cluster;
      int jcluster = box->Particles[it->j].cluster;
      if (isInner == 0 && icluster == jcluster) continue;
      if (isInner == 1 && icluster != jcluster) continue;

      // Since it is not possible to modify an element in a std::set
      // we break this restriction by defining and using the following pointer
      BreakableInterface* I = const_cast<BreakableInterface*>(std::addressof(*it));
      
      int igroup = box->Particles[it->i].group;
      int jgroup = box->Particles[it->j].group;
      
      if (isInner) {       
        I->kn = box->dataTable.get(box->idKnInnerBond, igroup, jgroup);
        I->kt = box->dataTable.get(box->idKtInnerBond, igroup, jgroup);
        I->fn0 = box->dataTable.get(box->idFn0InnerBond, igroup, jgroup);
        I->ft0 = box->dataTable.get(box->idFt0InnerBond, igroup, jgroup);
        I->power = box->dataTable.get(box->idPowInnerBond, igroup, jgroup);
      } else {
        I->kn = box->dataTable.get(box->idKnOuterBond, igroup, jgroup);
        I->kt = box->dataTable.get(box->idKtOuterBond, igroup, jgroup);
        I->kr = box->dataTable.get(box->idKrOuterBond, igroup, jgroup);
        I->fn0 = box->dataTable.get(box->idFn0OuterBond, igroup, jgroup);
        I->ft0 = box->dataTable.get(box->idFt0OuterBond, igroup, jgroup);
        I->mom0 = box->dataTable.get(box->idMom0OuterBond, igroup, jgroup);
        I->power = box->dataTable.get(box->idPowOuterBond, igroup, jgroup);
      }
    }
  }
  box->paramsInInterfaces = 1;  // so that they are stored in the conf files (within interfaces)
}
