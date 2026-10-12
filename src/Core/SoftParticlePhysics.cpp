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

#include "SoftParticlePhysics.hpp"

SoftParticlePhysics::SoftParticlePhysics() {
  F = mat9r::unit();
  stress.reset();
  strain.reset();
  stressSum.reset();
}

void SoftParticlePhysics::setElasticity(double Young, double Poisson) {
  young = Young;
  poisson = Poisson;
  Cinv.set(Young, Poisson);
  C.set(Young, Poisson);
}

/**
 * @brief Makes the material viscoelastic (Kelvin-Voigt) with the given shear and bulk viscosities.
 *        Zero viscosities give back the elastic material.
 */
void SoftParticlePhysics::setViscosity(double ShearVis, double BulkVis) {
  shearViscosity = ShearVis;
  bulkViscosity = BulkVis;
  useVisco = (ShearVis > 0.0 && BulkVis > 0.0);
  if (useVisco) {
    Dinv.set(ShearVis, BulkVis);
  }
}

/**
 * @brief Makes the material viscoelastic with a viscosity proportional to the stiffness, D = tau C.
 *        The elasticity has to be set before.
 */
void SoftParticlePhysics::setViscousTime(double tau) {
  double G = young / (2.0 * (1.0 + poisson));
  double K = young / (3.0 * (1.0 - 2.0 * poisson));
  setViscosity(tau * G, tau * K);
}

/**
 * @brief Turns the accumulated stressSum into the mean stress of the body, and updates the deformation gradient.
 *
 * Elastic material: the stress increment since the last update gives a strain increment, so that the update
 * does not depend on the time step.
 * Viscoelastic material (Kelvin-Voigt): the strain is integrated explicitly, eps += dt Dinv (sigma - C eps).
 * In both cases, the strain increment is applied exactly as F <- exp(dEps) F.
 *
 * @param Q                 Orientation of the body
 * @param undeformedVolume  Volume of the body before its deformation (homothety included)
 * @param dt                Time increment since the last update (used by the viscoelastic material only)
 */
void SoftParticlePhysics::updateTransformation(const quat& Q, double undeformedVolume, double dt) {
  mat9r sig = stressSum / (undeformedVolume * F.det());
  sig.symmetrize();

  // Express the stress in the body frame
  mat9r R;
  Q.get_rot_matrix(R);
  mat9r stressBody = R.transposed() * sig * R;

  mat9r dEps;
  if (useVisco) {
    dEps = dt * Dinv.getStrainRate(stressBody - C.getStress(strain));
    strain += dEps;
  } else {
    dEps = Cinv.getStrain(stressBody - stress);
    strain = Cinv.getStrain(stressBody);
  }
  F = dEps.exp() * F;
  stress = stressBody;
}
