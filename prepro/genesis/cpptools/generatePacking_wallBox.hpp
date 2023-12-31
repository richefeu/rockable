#ifndef GENERATEPACKING_WALLBOX_HPP
#define GENERATEPACKING_WALLBOX_HPP

#include <iostream>

// Lengths LX, LY and LZ are inside (walls are placed outside this cube)
// The rank-order is xmin, xmax, ymin, ymax, zmin, zmax
// TODO: add an origin
int generatePacking_wallBox(std::ostream& os, int group, double LX, double LY, double LZ, double Rw) {
  using namespace std;

  // normal X
  os << "x-wall " << group << " 0 1 " << -Rw << " " << 0.5 * LY << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;
  os << "x-wall " << group << " 0 1 " << LX + Rw << " " << 0.5 * LY << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;

  // normal Y
  os << "y-wall " << group << " 0 1 " << 0.5 * LX << " " << -Rw << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;
  os << "y-wall " << group << " 0 1 " << 0.5 * LX << " " << LY + Rw << " " << 0.5 * LZ
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;

  // normal Z
  os << "z-wall " << group << " 0 1 " << 0.5 * LX << " " << 0.5 * LY << " " << -Rw
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;
  os << "z-wall " << group << " 0 1 " << 0.5 * LX << " " << 0.5 * LY << " " << LZ + Rw
     << "  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0 " << endl;

  return 6;
}

#endif /* end of include guard: GENERATEPACKING_WALLBOX_HPP */
