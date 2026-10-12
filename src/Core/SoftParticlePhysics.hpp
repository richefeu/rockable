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

#ifndef SOFTPARTICLEPHYSICS_HPP
#define SOFTPARTICLEPHYSICS_HPP

#include <cstddef>
#include <memory>

// toofus
#include "Compliance.hpp"
#include "Rigidity.hpp"
#include "mat9.hpp"
#include "quat.hpp"

#include "VisDamping.hpp"

/// The homogeneous deformation of a soft particle.
///
/// The deformation gradient F maps the shape (body frame) onto the deformed body, so that a vertex is placed at
/// pos + Q * (F * (homothety * vertex)). The Minkowski radius is not deformed.
/// The stress and the strain are expressed in the body frame.
///
/// The material is linear elastic, or viscoelastic (Kelvin-Voigt: sigma = C eps + D deps/dt) when a viscosity is set.
class SoftParticlePhysics {
 public:
  double young{0.0};          ///< Young modulus
  double poisson{0.0};        ///< Poisson ratio
  double shearViscosity{0.0}; ///< Shear viscosity (0 for an elastic material)
  double bulkViscosity{0.0};  ///< Bulk viscosity (0 for an elastic material)
  bool useVisco{false};       ///< Viscoelastic (Kelvin-Voigt) update instead of the elastic one

  Compliance Cinv;  ///< Elastic compliance built from young and poisson
  Rigidity C;       ///< Elastic stiffness built from young and poisson
  VisDamping Dinv;  ///< Inverse viscosity built from shearViscosity and bulkViscosity

  mat9r F;          ///< Deformation gradient (body frame)
  mat9r stress;     ///< Mean stress in the body (body frame), at the last update
  mat9r strain;     ///< Strain of the body (body frame), at the last update
  mat9r stressSum;  ///< Sum of the f x branch products (global frame), accumulated during the force computation

  SoftParticlePhysics();  // Ctor

  void setElasticity(double Young, double Poisson);
  void setViscosity(double ShearVis, double BulkVis);
  void setViscousTime(double tau);
  void updateTransformation(const quat& Q, double undeformedVolume, double dt);
};

/// Owning pointer to the SoftParticlePhysics of a particle.
/// Copying a particle copies its deformation state, so that each particle has its own one.
class SoftParticlePhysicsPtr {
  std::unique_ptr<SoftParticlePhysics> ptr;

 public:
  SoftParticlePhysicsPtr() = default;
  SoftParticlePhysicsPtr(const SoftParticlePhysicsPtr& other)
      : ptr(other.ptr ? std::make_unique<SoftParticlePhysics>(*other.ptr) : nullptr) {}
  SoftParticlePhysicsPtr(SoftParticlePhysicsPtr&&) = default;
  SoftParticlePhysicsPtr& operator=(const SoftParticlePhysicsPtr& other) {
    if (this != &other) {
      ptr = other.ptr ? std::make_unique<SoftParticlePhysics>(*other.ptr) : nullptr;
    }
    return *this;
  }
  SoftParticlePhysicsPtr& operator=(SoftParticlePhysicsPtr&&) = default;

  void create() { ptr = std::make_unique<SoftParticlePhysics>(); }
  SoftParticlePhysics* operator->() const { return ptr.get(); }
  SoftParticlePhysics& operator*() const { return *ptr; }
  bool operator==(std::nullptr_t) const { return ptr == nullptr; }
  bool operator!=(std::nullptr_t) const { return ptr != nullptr; }
};

#endif /* end of include guard: SOFTPARTICLEPHYSICS_HPP */
