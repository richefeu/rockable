#!/bin/sh
# Gathers the results of the sweep: for each computation listed in
# all_inputs.txt, prints en2 and mu (read in its input file), the final
# position x of the grain of rice (last conf-file) and the kinetic energy
# (translation + rotation) at the end of the computation.
# Usage: sh results.sh [list of the input files, default all_inputs.txt]
export LC_ALL=C
echo "# en2 mu x_rice kinetic_energy"
while read -r input; do
  dir=$(dirname "$input")
  last=$(ls "$dir" | sed -n 's/^conf\([0-9]*\)$/\1/p' | sort -n | tail -1)
  [ -n "$last" ] || { echo "# $dir: no conf-file" >&2; continue; }
  ek=$(tail -1 "$dir/kineticEnergy.txt" | awk '{print $2 + $3}')
  awk -v ek="$ek" '
    $1 == "en2Contact" && $2 == 0 && $3 == 1 { en2 = $4 }
    $1 == "muContact" && $2 == 0 && $3 == 1 { mu = $4 }
    $1 == "Particles" { n = $2; getline; for (i = 0; i < n; i++) { getline; x[i] = $5 } }
    END { printf "%s %s %.3f %.3g\n", en2, mu, x[1], ek }' "$dir/conf$last"
done < "${1:-all_inputs.txt}"
