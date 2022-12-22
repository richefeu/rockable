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

#include "PreproCommand_setVariableStickParams.hpp"
#include "Core/Rockable.hpp"

setVariableStickParams::setVariableStickParams() { }

void setVariableStickParams::addCommand() {  
  box->parser.kwMap["setVariableStickParams"] = [this](std::istream& conf) {
    int timeSeededInt = 0;
    conf >> this->paramName >> this->isInnerStr >> this->lambda >> this->m >> timeSeededInt;
    this->timeSeeded = (bool)timeSeededInt;
    exec();
  };
}

/**
   Usage (in an input file):
   variableStickParams paramName inner/outer lambdaValue mValue 0/1
   example: variableStickParams fn0 inner 80 9 1
*/
void setVariableStickParams::exec() {
  std::default_random_engine generator;
  if (timeSeeded == true) {
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }
  std::weibull_distribution<> distribution(m, lambda);

  int shift = 0;
  if (paramName == "kn")
    shift = 0;
  else if (paramName == "kt")
    shift = 1;
  else if (paramName == "kr")
    shift = 2;
  else if (paramName == "fn0")
    shift = 3;
  else if (paramName == "ft0")
    shift = 4;
  else if (paramName == "mom0")
    shift = 5;

  int isInner = 0;
  if (isInnerStr == "inner") isInner = 1;

  for (size_t i = 0; i < box->Interfaces.size(); i++) {
    for (auto it = box->Interfaces[i].begin(); it != box->Interfaces[i].end(); ++it) {
      int icluster = box->Particles[it->i].cluster;
      int jcluster = box->Particles[it->j].cluster;
      if (isInner == 0 && icluster == jcluster) continue;
      if (isInner == 1 && icluster != jcluster) continue;

      double* seek = (double*)&(it->kn);
      *(seek + shift) = distribution(generator);
    }
  }

  box->paramsInInterfaces = 1;  // so that they are stored in the conf files (within interfaces)
}
