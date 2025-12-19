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
//
#ifndef PERIODICCELL_HPP
#define PERIODICCELL_HPP

#include "mat9.hpp"

class PeriodicCell {
 public:
  mat9r h;      ///< Matrix that hold the cell geometry (each column is a vector that defines a side)
  mat9r hinv;   ///< inverse of matrix h (recomputed each time h is updated)
  mat9r vh;     ///< Velocities of the collective DoF
  mat9r ah;     ///< Acceleration of the collective DoF
  mat9r Sig;    ///< Mean stress matrix over the cell
  double mass;  ///< Mass of the periodic cell (a somehow fictive data)

  PeriodicCell();
  PeriodicCell(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY, double ZZ);
  void Define(double XX, double XY, double XZ, double YX, double YY, double YZ, double ZX, double ZY, double ZZ);
  void precomputeInverse();
  vec3r getBranchCorrection(const vec3r& ipos, const vec3r& jpos);
  void forceToStayInside(vec3r& s);
};

#endif /* end of include guard: PERIODICCELL_HPP */
