// Combined real-world example

#include "simuFileManips.hpp"
#include <vector>

int main() {
  std::vector<double> dtValues{1e-4, 1e-3};

  for (double dt : dtValues) {
    SFManip()
        .createFolder("example7/run_%g", dt)
        .setOutputFilename("config.txt")
        .silentIfNotFound()

        .read("input.txt")

        // core parameters
        .replaceAllStartingWith("dt ", "dt %g", dt)
        .replaceAllStartingWith("steps ", "steps %d", 10000)

        // tweak behavior
        .insertAfter("solver ", "tolerance 1e-8")

        // toggle debug if present
        .findLineStartingWith("debug ")
        .replaceInLine("false", "true")

        // metadata
        .appendLine("# generated automatically")

        .saveInFolder();
  }
}