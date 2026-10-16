#!/bin/sh
# Energy balance of the shot: the kinetic energy of the ball before the impact,
# the energy needed to break all the interfaces, the energy spent in the cracks,
# and the kinetic energy of the ball and of the cells at the end.
# An interface breaks when its elastic energy exceeds 2 A Gc (A its area, Gc the
# fracture energy, 2 for the two faces of the crack).
# usage: sh energies.sh [conf-file]          (conf40 by default)
export LC_ALL=C
CONF=${1:-conf40}

# kinetic energies of the ball (group 2) and of the cells (group 1), and the
# energy 2 A Gc summed over the interfaces of a conf-file; the volume of each
# shape is read in shapes.txt
energy() {
  awk '
    FNR == NR { if ($1 == "name") s = $2; if ($1 == "volume") V[s] = $2; next }
    $1 == "gcInnerBond" && $2 == 1 && $3 == 1 { Gc = $4 }
    $1 == "density" { rho[$2] = $3 }
    /^Particles/ { n = $2; reading = 1; next }
    reading && /^#/ { next }
    reading && n > 0 {
      n--
      e = 0.5 * rho[$2] * V[$1] * $4 * $4 * $4 * ($8 * $8 + $9 * $9 + $10 * $10)
      if ($2 == 2) Eb += e
      if ($2 == 1) Ec += e
      if (n == 0) reading = 0
      next
    }
    $1 == "Interfaces" { m = $2; next }
    m > 0 { Ef += 2.0 * $5 * Gc; m-- }
    END { printf "%g %g %g\n", Eb, Ec, Ef }' shapes.txt "$1"
}

set -- $(energy conf0) $(energy "$CONF")
awk -v Eb0="$1" -v Ef0="$3" -v Eb="$4" -v Ec="$5" -v Ef="$6" -v conf="$CONF" 'BEGIN {
  printf "kinetic energy of the ball before the impact  %8.3f J\n", Eb0
  printf "energy to break all the interfaces            %8.3f J\n", Ef0
  printf "energy spent in the cracks (%s)           %8.3f J  (%.1f%% of the kinetic energy)\n", conf, Ef0 - Ef, 100 * (Ef0 - Ef) / Eb0
  printf "kinetic energy left in the ball (%s)      %8.3f J\n", conf, Eb
  printf "kinetic energy of the cells (%s)          %8.3f J\n", conf, Ec
}'
