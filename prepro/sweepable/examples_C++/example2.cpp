#include "simuFileManips.hpp"
#include <vector>

int main() {
  std::vector<double> dtValues{1e-5, 5e-5, 1e-4};

  // Editing multiple parameters in one file
  for (double dt : dtValues) {
    SFManip()
        .createFolder("example2")
        .read("input.txt")
        .replaceAllStartingWith("dt ", "dt %f", dt)
        .replaceAllStartingWith("viscosity ", "viscosity %f", 1e-3)
        .replaceAllStartingWith("steps ", "steps %d", 10000)
        .saveInFolder("example2/dt%g", dt);
  }
}