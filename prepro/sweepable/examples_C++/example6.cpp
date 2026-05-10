#include "simuFileManips.hpp"

int main() {

  // Silent mode (robust scripts
  SFManip()
      .silentIfNotFound()
      .read("input.txt")
      .replaceAllStartingWith("optional_param ", "optional_param %d", 42)
      .insertAfter("maybe_here ", "extra_line true")
      .saveInFolder("example6/robust_run");
}