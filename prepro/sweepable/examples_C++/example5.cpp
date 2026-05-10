#include "simuFileManips.hpp"

int main() {

  // Custom output filename
  SFManip()
      .createFolder("example5/runA")
      .setOutputFilename("config.txt")
      .read("input.txt")
      .replaceAllStartingWith("dt ", "dt %f", 0.01)
      .saveInFolder();
}