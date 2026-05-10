#include "simuFileManips.hpp"
#include <vector>

int main() {
  std::vector<double> dtValues{1e-5, 5e-5, 1e-4};
  std::vector<int> maxIterValues{200, 400};

  SFManip::clearCollection();

  // Parameter sweep
  for (double dt : dtValues) {
    for (int maxIter : maxIterValues) {
      
      SFManip()
          .createFolder("example1/dt_%g/maxIter_%d", dt, maxIter)
          .read("input.txt")
          .replaceAllStartingWith("dt ", "dt %g", dt)
          .replaceAllStartingWith("max_iter ", "max_iter %d", maxIter)
          .appendLine("# auto-generated")
          .saveInFolder();
    }
  }
  
  SFManip::saveCollection("example1/lst-input.txt");
}

// Running commands in parallel from a file using xargs:
//
// Goal:
//   For each line in a file, run: ./compute <line>
//   with at most N processes running in parallel.
//
// Portable solution (works on Linux + macOS/BSD):
//   xargs -P N -n 1 ./run < lst-commands.txt
//
// Example with echo:
//   xargs -P 3 -n 1 echo "../../BUILD/rockable" < lst-commands.txt
//
// Notes:
// - "-P N" controls the number of parallel jobs.
// - "-n 1" ensures one line = one command.
// - Input is read via "< lst-commands.txt" (since "-a" is GNU-only).
// - Output order is not guaranteed due to parallel execution.

