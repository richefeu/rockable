#!/bin/sh
# Solid fraction of the sample of isoComprSpheres, during the compression.
#
# It reads the positions of the three moving walls in xmax.txt, ymax.txt and
# zmax.txt (written by the data extractors, see dataExtractors.txt), and the
# volume of the spheres in a conf-file (conf0 by default).
#
#   sh solidFraction.sh [conf-file]
#
# Columns of the output: time, Lx, Ly, Lz (inner size of the box), solid fraction

# Numbers are read and written with a decimal point, whatever the language settings
export LC_ALL=C

CONF=${1:-conf0}

# Volume of the spheres: the shape 'Sphere' has a unit radius, so the radius of
# each sphere is its homothety (4th column of its line in the conf-file)
VS=$(awk '$1 == "Sphere" { V += 4.0 / 3.0 * 3.14159265358979 * $4 * $4 * $4 } END { print V }' "$CONF")

# Each line of a TrackBody file has 20 columns: the time, then the position.
# Side by side, the x position of xmax is column 2, the y position of ymax is
# column 23 and the z position of zmax is column 44. The inner face of the fixed
# walls is at 0, and the moving walls have a radius of 0.001.
paste xmax.txt ymax.txt zmax.txt | awk -v Vs="$VS" '
  BEGIN { print "# time  Lx  Ly  Lz  solid_fraction" }
  NF >= 44 {
    Lx = $2 - 0.001; Ly = $23 - 0.001; Lz = $44 - 0.001
    printf "%.4f %.5f %.5f %.5f %.4f\n", $1, Lx, Ly, Lz, Vs / (Lx * Ly * Lz)
  }'
