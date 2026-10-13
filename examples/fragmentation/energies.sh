#!/bin/sh
# Energy balance of the impact: the kinetic energy of the block before the
# impact, and the energy spent in breaking interfaces.
#
#   sh energies.sh [conf-file]          (conf20 by default)
#
# An interface breaks when its elastic energy exceeds 2 A Gc (A its area, Gc
# the fracture energy, 2 for the two faces of the crack). The energy needed to
# break all of them is the sum of 2 A Gc over the interfaces of conf0, and the
# energy spent is the same sum over the interfaces that are missing in the
# conf-file given.

# Numbers are read and written with a decimal point, whatever the language settings
export LC_ALL=C

CONF=${1:-conf20}

# Volume of the shape Cell, in shapes.txt
VCELL=$(awk '/^name / { n = $2 } /^volume / && n == "Cell" { print $2 }' shapes.txt)

# In a conf-file: the line 'gcInnerBond 1 1 <Gc>', the lines of the cells
# (name, group, cluster, homothety, position, then the velocity in columns 8-10),
# and the block 'Interfaces' (the area is the 5th column of each line)
energy() {
  awk -v Vc="$VCELL" '
    $1 == "gcInnerBond" && $2 == 1 && $3 == 1 { Gc = $4 }
    $1 == "density" && $2 == 1 { rho = $3 }
    $1 == "Cell" { Ec += 0.5 * rho * Vc * $4 * $4 * $4 * ($8 * $8 + $9 * $9 + $10 * $10) }
    $1 == "Interfaces" { n = $2; next }
    n > 0 { Ef += 2.0 * $5 * Gc; n-- }
    END { printf "%g %g\n", Ec, Ef }' "$1"
}

set -- $(energy conf0) $(energy "$CONF")
awk -v Ec0="$1" -v Ef0="$2" -v Ec="$3" -v Ef="$4" -v conf="$CONF" 'BEGIN {
  printf "kinetic energy before the impact      %8.3f J\n", Ec0
  printf "energy to break all the interfaces    %8.3f J\n", Ef0
  printf "energy spent in the cracks (%s)   %8.3f J  (%.1f%% of the kinetic energy)\n", conf, Ef0 - Ef, 100 * (Ef0 - Ef) / Ec0
  printf "kinetic energy left (%s)          %8.3f J\n", conf, Ec
}'
