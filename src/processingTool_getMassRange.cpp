#include "processingTool_getMassRange.hpp"

/**
   @brief  Get the range of masses

   @param[out]  massMin  The lightest mass (for particles id in range[first last])
   @param[out]  massMax  The heaviest mass (for particles id in range[first last])
   @param[in]   first    Smallest ID of particles (default value is 0)
   @param[in]   last     Largest ID of particles (default value corresponds to the last particle)
*/
void getMassRange(Rockable *box, double& massMin, double& massMax, size_t first, size_t last) {
  if (last == 0) last = box->Particles.size() - 1;
  massMin = massMax = box->Particles[first].mass;
  for (size_t i = first + 1; i <= last; i++) {
    if (box->Particles[i].mass < massMin) massMin = box->Particles[i].mass;
    if (box->Particles[i].mass > massMax) massMax = box->Particles[i].mass;
  }
}
