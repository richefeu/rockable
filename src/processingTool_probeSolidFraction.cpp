#include "processingTool_probeSolidFraction.hpp"

/**
 * @brief xxx
 *
 * @param aabb    The probe as an Axis Aligned Bounding Box
 * @param MCnstep Number of Monte Carlo steps
 * @return double Estimate of the solid fraction inside the probe
 */
double probeSolidFraction(Rockable *box, AABB& aabb, size_t MCnstep) {
  if (MCnstep == 0) return -1.0;

  // select the concerned particles
  OBB zone;
  zone.center = 0.5 * (aabb.min + aabb.max);
  zone.extent.set(0.5 * (aabb.max.x - aabb.min.x), 0.5 * (aabb.max.y - aabb.min.y), 0.5 * (aabb.max.z - aabb.min.z));
  std::vector<size_t> pid;
  for (size_t i = 0; i < box->Particles.size(); ++i) {
    box->Particles[i].updateObb();
    if (zone.intersect(box->Particles[i].obb)) {
      pid.push_back(i);
    }
  }

  vec3r pt3;
  std::vector<double> vv(3);
  Mth::sobolSequence(-3, vv);  // Initialize the Sobol sequence
  size_t count = 0;
  for (size_t imc = 0; imc < MCnstep; ++imc) {
    Mth::sobolSequence(3, vv);
    pt3.set(aabb.min.x + vv[0] * (aabb.max.x - aabb.min.x), aabb.min.y + vv[1] * (aabb.max.y - aabb.min.y),
            aabb.min.z + vv[2] * (aabb.max.z - aabb.min.z));

    bool inSolid = false;
    for (size_t ii = 0; ii < pid.size(); ii++) {
      size_t i = pid[ii];
      vec3r ptTest = pt3 - box->Particles[i].pos;
      quat Qinv = box->Particles[i].Q.get_conjugated();
      ptTest = Qinv * ptTest;
      ptTest /= box->Particles[i].homothety;

      if (box->Particles[i].shape->inside(ptTest)) {
        inSolid = true;
        break;
      }
    }
    if (inSolid) count++;
  }

  return ((double)count / (double)MCnstep);
}