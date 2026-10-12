#!/bin/sh
# Solid fraction of the periodic sample, read in the conf-files.
# Prints: t  Lx  Ly  Lz  solid-fraction
# The cell volume is det(h); the spheres have a unit-radius shape, so the
# radius of each one is its homothety.
# Usage: sh solidFraction.sh [first last step]   (default 0 100 1)
export LC_ALL=C
first=${1:-0}; last=${2:-100}; step=${3:-1}
i=$first
while [ $i -le $last ] && [ -f conf$i ]; do
  awk '
    $1 == "t" { t = $2 }
    $1 == "h" { a=$2; b=$3; c=$4; d=$5; e=$6; f=$7; g=$8; hh=$9; k=$10
                V = a*(e*k - f*hh) - b*(d*k - f*g) + c*(d*hh - e*g)
                Lx = a; Ly = e; Lz = k }
    $1 == "Particles" { n = $2; getline; for (p = 0; p < n; p++) { getline; vs += 4.0 / 3.0 * 3.141592653589793 * $4 ^ 3 } }
    END { printf "%g %.5f %.5f %.5f %.4f\n", t, Lx, Ly, Lz, vs / V }' conf$i
  i=$((i + step))
done
