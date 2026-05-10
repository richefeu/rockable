#include "simuFileManips.hpp"

int main() {

  // Insert configuration blocks
  SFManip()
      .read("input.txt")
      .insertAfter("solver ", "tolerance 1e-6")
      .insertAfter("solver ", "max_iter 500")
      .saveInFolder("example4/with_solver_settings");
}