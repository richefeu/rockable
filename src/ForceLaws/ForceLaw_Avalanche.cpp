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
#include "ForceLaw_Avalanche.hpp"

Avalanche::Avalanche() { }

void Avalanche::init() {
  box->idKnContact = box->dataTable.add("knContact");
  box->idEn2Contact = box->dataTable.add("en2Contact");
  box->idKtContact = box->dataTable.add("ktContact");
  box->idMuContact = box->dataTable.add("muContact");
  box->idKrContact = box->dataTable.add("krContact");
  box->idMurContact = box->dataTable.add("murContact");
}

/**
    @brief   Force-law used for rock avalanches at Laboratoire 3SR
    @return  Return true if the interaction is active (ie. with a non-zero force)
*/
bool Avalanche::computeInteraction(Interaction& I) {
  if (I.dn > 0.0) {
    I.fn = 0.0;
    I.ft.reset();
    I.mom.reset();
    return false;
  }

  // getInteractingGroups(I, g1, g2);
  int g1 = box->Particles[I.i].group;
  int g2 = box->Particles[I.j].group;
  double kn = box->dataTable.get(box->idKnContact, g1, g2);
  double en2 = box->dataTable.get(box->idEn2Contact, g1, g2);
  double kt = box->dataTable.get(box->idKtContact, g1, g2);
  double mu = box->dataTable.get(box->idMuContact, g1, g2);
  double kr = box->dataTable.get(box->idKrContact, g1, g2);
  double mur = box->dataTable.get(box->idMurContact, g1, g2);

  if (box->ctcPartnership.getWeight != nullptr) {
    double w = box->ctcPartnership.getWeight(I);
    kn *= w;
    kt *= w;
    kr *= w;
  }

  // === Normal force
  if (I.prev_dn > 0.0) I.prev_dn = 0.0;
  double delta_Dn = I.dn - I.prev_dn;
  if (delta_Dn > 0.0)
    I.fn = -kn * en2 * I.dn;  // Unloading
  else if (delta_Dn < 0.0)
    I.fn += -kn * delta_Dn;  // Loading
  if (I.fn < 0.0) I.fn = 0.0;

  // === Tangential force (friction)
  vec3r vt = I.vel - (I.vel * I.n) * I.n;
#ifdef FT_CORR
  vec3r ft_corr = I.ft;
  ft_corr -= cross(ft_corr, cross(I.prev_n, I.n));
  ft_corr -= cross(ft_corr, (box->dt_2 * (box->Particles[I.i].vrot + box->Particles[I.j].vrot) * I.n) * I.n);
  I.ft = ft_corr + kt * (vt * box->dt);
#else
  I.ft += kt * (box->dt * vt);
#endif
  double threshold_ft = fabs(mu * I.fn);  // even without fabs the value should be positive
  double ft_square = I.ft * I.ft;
  if (ft_square > 0.0 && ft_square > threshold_ft * threshold_ft) I.ft *= threshold_ft / sqrt(ft_square);
  // Remark: in fact, the test (ft * ft > 0.0) means that ft_square is NOT null,
  // because ft * ft >= 0 by definition.

  // === Resistant moment
  I.mom += kr * (box->Particles[I.j].vrot - box->Particles[I.i].vrot) * box->dt;
  vec3r branch;
  if (I.i < box->nDriven) {  // j is the free body (a rock block)
    branch = I.pos - box->Particles[I.j].pos;
    double r = (branch * box->Particles[I.j].vrot) / (box->Particles[I.j].vrot * box->Particles[I.j].vrot);
    branch -= r * box->Particles[I.j].vrot;
  } else {  // i is the free body (a rock block)
    branch = I.pos - box->Particles[I.i].pos;
    double r = (branch * box->Particles[I.i].vrot) / (box->Particles[I.i].vrot * box->Particles[I.i].vrot);
    branch -= r * box->Particles[I.i].vrot;
  }
  double threshold_mom = fabs(mur * norm(branch) * I.fn);  // even without fabs, the value should
                                                           // be positive
  double mom_square = I.mom * I.mom;
  if (mom_square > 0.0 && mom_square > threshold_mom * threshold_mom) I.mom *= threshold_mom / sqrt(mom_square);

  return true;
}
