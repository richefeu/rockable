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
#include "ForceLaw_geoVisc.hpp"

GeoVisc::GeoVisc() {}

void GeoVisc::init() {
  box->idKnContact = box->dataTable.add("knContact");
  box->idEn2Contact = box->dataTable.add("en2Contact");
  box->idKtContact = box->dataTable.add("ktContact");
  box->idMuContact = box->dataTable.add("muContact");
  box->idKrContact = box->dataTable.add("krContact");
  box->idMurContact = box->dataTable.add("murContact");
}

/**
    @brief   Based on the default law, this law is designed for slow landslides.
             It includes Linear normal repulsion, normal viscosity, 'Coulomb-like' viscuous Friction, 
             and moment resistance. The specific feature is to make the friction force viscuous while still using 
             the Coulomb cone yielding.
    @return  Return true if the interaction is active (ie. with a non-zero force)
*/
bool GeoVisc::computeInteraction(Interaction& I) {
  if (I.dn > 0.0) {
    I.fn = 0.0;
    I.ft.reset();
    I.mom.reset();
    return false;
  }

  int g1 = box->Particles[I.i].group;
  int g2 = box->Particles[I.j].group;
  double kn = box->dataTable.get(box->idKnContact, g1, g2);
  double viscT = box->dataTable.get(box->idViscTContact, g1, g2);
  double mu = box->dataTable.get(box->idMuContact, g1, g2);
  double kr = box->dataTable.get(box->idKrContact, g1, g2);
  double mur = box->dataTable.get(box->idMurContact, g1, g2);
  double damp = I.damp;  // It has been precomputed

  if (box->ctcPartnership.getWeight != nullptr) {
    double w = box->ctcPartnership.getWeight(I);
    kn *= w;
    //kt *= w;
    kr *= w;
    damp *= sqrt(w);
  }

  // === Normal force (elatic contact + viscous damping)
  double vn = I.vel * I.n;
  double fne = -kn * I.dn;

  /*
  if (box->preventCrossingLength > 0.0) {
    fne *= box->preventCrossingLength / (box->preventCrossingLength + I.dn);
  }
  */

  double fnv = damp * vn;
  I.fn = fne + fnv;
  if (I.fn < 0.0) {
    I.fn = 0.0;  // Because viscous damping can make the normal force negative
  }

  // === Tangential force (viscuous friction)
  vec3r vt = I.vel - vn * I.n;
  I.ft = viscT * vt;

#ifdef ALLOW_NEGATIVE_THRESHOLDS
  double threshold = fabs(mu * I.fn);
#else
  double threshold = fabs(mu * fne);
#endif

  double ft_square = I.ft * I.ft;
  if (ft_square > 0.0 && ft_square >= threshold * threshold) {
    I.ft *= threshold / sqrt(ft_square);
    // Remark: in fact, the test (ft_square > 0.0) means that ft_square is not
    // zero, because ft_square >= 0 by definition.
  }

  // === Resistant moment
  if (kr > 0.0) {
    I.mom += kr * (box->Particles[I.j].vrot - box->Particles[I.i].vrot) * box->dt;
    double threshold_mom = fabs(mur * I.fn);  // in this case mur is a *length*
    double mom_square = I.mom * I.mom;
    if (mom_square > 0.0 && mom_square >= threshold_mom * threshold_mom) I.mom *= threshold_mom / sqrt(mom_square);
  }

  return true;
}
