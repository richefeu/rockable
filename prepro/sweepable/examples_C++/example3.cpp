#include "simuFileManips.hpp"

int main() {
  
  // Modify part of a line (not full replacement)
  SFManip()
      .read("input.txt")
      .findLineStartingWith("output ")
      .replaceInLine("false", "true")
      .saveInFolder("example3/debug_run");
}