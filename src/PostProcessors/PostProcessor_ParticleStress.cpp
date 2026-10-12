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

#include "PostProcessor_ParticleStress.hpp"
#include "Core/Rockable.hpp"

// Registered in Rockable::ExplicitRegistrations

ParticleStress::ParticleStress() { }

void ParticleStress::read(std::istream& is) {
  kwParser parser;
  parser.kwMap["Volume"] = __GET__(istr, Volume);
  parser.kwMap["ConfVolumes"] = __DO__(istr) {
    size_t nb;
    istr >> nb;
    for (size_t i = 0 ; i < nb ; i++) {
      int iconf;
      double v;
      istr >> iconf >> v;
      ConfVolumes[iconf] = v;
    }
  };
	
  parser.parse(is);
  
  // example:
  // firstConf 0 
  // lastConf 2
  // stepConf 1
  // PostProcessor ParticleStress
  // Volume 1.0 # the default total volume
  // ConfVolumes 3 # volume corresponding to each conf (iconf value in the file)
  // 0 0.0012
  // 1 0.00121
  // 2 0.00123
  
}

void ParticleStress::init() {

}

void ParticleStress::end() {

}

void ParticleStress::exec() {
  char fname[256];
  sprintf(fname, "ParticleStress-conf%d.txt", box->iconf);
  std::ofstream file(fname);

  // Vector from the contact point to the centre of the particle (nearest image with a periodic cell)
  auto lever = [this](const vec3r& xp, const vec3r& xc) {
    vec3r d = xp - xc;
#ifdef ROCKABLE_ENABLE_PERIODIC
    if (box->usePeriodicCell == 1) {
      vec3r s = box->Cell.hinv * d;
      s.x -= std::round(s.x);
      s.y -= std::round(s.y);
      s.z -= std::round(s.z);
      d = box->Cell.h * s;
    }
#endif
    return d;
  };

  std::vector<mat9r> M(box->Particles.size());
  for (size_t ibody = 0; ibody < box->Particles.size(); ibody++) {
    for (auto it = box->Interactions[ibody].begin(); it != box->Interactions[ibody].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));

      size_t i = I->i;
      size_t j = I->j;

      // The normal n goes from j to i, so that fn * n + ft is the force received by i.
      // M_p = sum over the contacts c of F_c (x) (x_p - x_c), with F_c the force received by p: it is positive
      // in compression, and it does not depend on the origin, even when the particle is not in equilibrium
      // (with the absolute position x_c, the moment of the resultant force was added to it).
      vec3r f = I->fn * I->n + I->ft;

      if (i >= box->nDriven) {  // the tensorial moment of walls is not computed
        M[i] += dyadic_product(f, lever(box->Particles[i].pos, I->pos));
      }

      if (j >= box->nDriven) {  // the tensorial moment of walls is not computed
        M[j] += dyadic_product(-f, lever(box->Particles[j].pos, I->pos));
      }
    }
  }

  mat9r Mtotal;
  double Vtot = 1.0; 
  if (ConfVolumes.find(box->iconf) != ConfVolumes.end())
    Vtot = ConfVolumes[box->iconf];
  else
    Vtot = Volume;
  __SHOW(box->iconf);
  __SHOW(Vtot);
  
  // This can be read in the postpro command file.
  // But in case of multiple conf-files, we need many
  vec3r bodyAxis(0.0, 0.0, 1.0);
  
  for (size_t ibody = 0; ibody < box->Particles.size(); ibody++) {
    vec3r u = box->Particles[ibody].Q * bodyAxis;
    
    // Two unit vectors perpendicular to u, built from an axis that is not parallel to u
    vec3r helper = (fabs(u.y) < 0.9) ? vec3r(0.0, 1.0, 0.0) : vec3r(1.0, 0.0, 0.0);
    vec3r ut1 = cross(u, helper);
    ut1.normalize();
    vec3r ut2 = cross(ut1, u);
    double a = (M[ibody] * u) * u;
    double r1 = (M[ibody] * ut1) * ut1;
    double r2 = (M[ibody] * ut2) * ut2;
    double r = sqrt(r1*r1+r2*r2);
    
    mat9r P = dyadic_product(u, u);
    mat9r Q = mat9r::unit() - P;
   
    mat9r Ma = P * M[ibody];
    mat9r Mr = Q * M[ibody];
      
    file << M[ibody] << ' '
      << a << ' ' << r << ' '
      << Ma.trace() << ' ' << Mr.trace() << std::endl;
    Mtotal += M[ibody];
  }
  
  std::cout << "total: " << Mtotal / Vtot << std::endl;
}
