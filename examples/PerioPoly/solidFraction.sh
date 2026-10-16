#!/bin/sh
# Solid fraction of the periodic sample, read in the conf-files.
# Prints: t  Lx  Ly  Lz  solid-fraction
# The volume of each particle is the volume of its shape (read in shape.shp)
# times its homothety to the cube; the volume of the cell is det(h).
# Usage: sh solidFraction.sh [first last step]   (default 0 30 1)
export LC_ALL=C
first=${1:-0}; last=${2:-30}; step=${3:-1}
i=$first
while [ $i -le $last ] && [ -f conf$i ]; do
  awk '
    FNR == NR { if ($1 == "name") name = $2; if ($1 == "volume") vol[name] = $2; next }
    $1 == "t" { t = $2 }
    $1 == "h" { V = $2*($6*$10 - $7*$9) - $3*($5*$10 - $7*$8) + $4*($5*$9 - $6*$8)
                Lx = $2; Ly = $6; Lz = $10 }
    $1 == "Particles" { n = $2; getline; for (p = 0; p < n; p++) { getline; vs += vol[$1] * $4 ^ 3 } }
    END { printf "%g %.4f %.4f %.4f %.4f\n", t, Lx, Ly, Lz, vs / V }' shape.shp conf$i
  i=$((i + step))
done
