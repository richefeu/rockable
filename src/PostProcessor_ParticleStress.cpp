// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#include "factory.hpp"
#include "kwParser.hpp"

#include "PostProcessor_ParticleStress.hpp"
#include "Rockable.hpp"

static Registrar<PostProcessor, ParticleStress> registrar("ParticleStress");

ParticleStress::ParticleStress() { }

void ParticleStress::read(std::istream& is) {
  kwParser parser;
  parser.kwMap["Volume"] = __GET__(is, Volume);
  parser.kwMap["ConfVolumes"] = __DO__(is) {
    size_t nb;
    is >> nb;
    for (size_t i = 0 ; i < nb ; i++) {
      int iconf;
      double v;
      is >> iconf >> v;
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

  std::vector<mat9r> M(box->Particles.size());
  for (size_t ibody = 0; ibody < box->Particles.size(); ibody++) {
    for (auto it = box->Interactions[ibody].begin(); it != box->Interactions[ibody].end(); ++it) {
      Interaction* I = const_cast<Interaction*>(std::addressof(*it));

      size_t i = I->i;
      size_t j = I->j;

      vec3r x = I->pos;
      vec3r f = -(I->fn * I->n + I->ft); // the sign - is for obtaining positive M for compression

      double sxx = f.x * x.x;
      double sxy = f.x * x.y;
      double sxz = f.x * x.z;
      double syx = f.y * x.x;
      double syy = f.y * x.y;
      double syz = f.y * x.z;
      double szx = f.z * x.x;
      double szy = f.z * x.y;
      double szz = f.z * x.z;

      if (i >= box->nDriven) { // the tensorial moment of walls is not computed
        M[i].xx += sxx;
        M[i].xy += sxy;
        M[i].xz += sxz;
        M[i].yx += syx;
        M[i].yy += syy;
        M[i].yz += syz;
        M[i].zx += szx;
        M[i].zy += szy;
        M[i].zz += szz;
      }

      if (j >= box->nDriven) { // the tensorial moment of walls is not computed
        M[j].xx -= sxx;
        M[j].xy -= sxy;
        M[j].xz -= sxz;
        M[j].yx -= syx;
        M[j].yy -= syy;
        M[j].yz -= syz;
        M[j].zx -= szx;
        M[j].zy -= szy;
        M[j].zz -= szz;
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
    
    vec3r y(0.0, 1.0, 0.0);
    vec3r ut1;
    if (fabs(u.x) < 1e-10 && fabs(u.z) < 1e-10) ut1 = cross(u, y); // u different from y
    else ut1.set(1.0, 0.0, 0.0);
    ut1.normalize();
    vec3r ut2 = cross (ut1,u);
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
