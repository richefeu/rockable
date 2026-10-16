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

#ifndef VISDAMPING_HPP
#define VISDAMPING_HPP

// toofus
#include "mat9.hpp"

/// Inverse of the viscosity tensor D of an isotropic Kelvin-Voigt material (written by Mukesh Singh Bisht).
///
/// It has the same form as the compliance:
///
///     Dinv1 Dinv2 Dinv2  0     0     0
///     Dinv2 Dinv1 Dinv2  0     0     0
///     Dinv2 Dinv2 Dinv1  0     0     0
///     0     0     0      Dinv3 0     0
///     0     0     0      0     Dinv3 0
///     0     0     0      0     0     Dinv3
///
/// With a shear viscosity tau*G and a bulk viscosity tau*K, D = tau*C (C being the elastic stiffness).
class VisDamping {
 public:
  double Dinv1{0.0}, Dinv2{0.0}, Dinv3{0.0};

  VisDamping() {}
  VisDamping(double ShearVis, double BulkVis) { set(ShearVis, BulkVis); }

  /**
   * @brief Set the inverse viscosity from the shear and bulk viscosities
   *
   * The strains are tensor components (not the engineering shear strains), so the shear term is 1/(2 ShearVis),
   * as 1/(2G) in the compliance.
   */
  void set(double ShearVis, double BulkVis) {
    Dinv1 = 1.0 / (9.0 * BulkVis) + 1.0 / (3.0 * ShearVis);
    Dinv2 = 1.0 / (9.0 * BulkVis) - 1.0 / (6.0 * ShearVis);
    Dinv3 = 1.0 / (2.0 * ShearVis);
  }

  /**
   * @brief Strain rate produced by a (viscous) stress
   */
  mat9r getStrainRate(const mat9r& stress) const {
    double Exx = Dinv1 * stress.xx + Dinv2 * (stress.yy + stress.zz);
    double Eyy = Dinv1 * stress.yy + Dinv2 * (stress.xx + stress.zz);
    double Ezz = Dinv1 * stress.zz + Dinv2 * (stress.xx + stress.yy);
    double Eyz = Dinv3 * stress.yz;
    double Exz = Dinv3 * stress.xz;
    double Exy = Dinv3 * stress.xy;
    // clang-format off
    return mat9r(
      Exx, Exy, Exz,
      Exy, Eyy, Eyz,
      Exz, Eyz, Ezz
    );
    // clang-format on
  }
};

#endif /* end of include guard: VISDAMPING_HPP */
