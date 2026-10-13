#!/bin/sh
# The vertical stretch of each soft bunny along the computation: the component
# F_yy of its deformation gradient, turned into the frame of the floor (1 at
# rest, 0.6 for a bunny squashed by 40%). One line per conf-file: the time, then
# F_yy for each soft particle, in the order of the section SoftParticles.
# The last line gives the smallest value of each column.
# usage: sh squash.sh > squash.txt
export LC_ALL=C
i=0
while [ -f conf$i ]; do
  awk '
  /^t / { t = $2 }
  /^Particles/ { n = $2; k = 0; reading = 1; next }
  reading && /^#/ { next }
  reading && n > 0 { w[k] = $14; x[k] = $15; y[k] = $16; z[k] = $17; k++; n--; if (n == 0) reading = 0; next }
  /^SoftParticles/ { m = $2; getline
    for (s = 0; s < m; s++) {
      getline
      id = $1
      for (c = 0; c < 9; c++) F[c] = $(6 + c)
      # the row y of the rotation matrix of the particle
      R[0] = 2 * (x[id] * y[id] + w[id] * z[id])
      R[1] = 1 - 2 * (x[id] * x[id] + z[id] * z[id])
      R[2] = 2 * (y[id] * z[id] - w[id] * x[id])
      f = 0
      for (a = 0; a < 3; a++) for (b = 0; b < 3; b++) f += R[a] * F[3 * a + b] * R[b]
      line = line sprintf(" %.4f", f)
    }
  }
  END { printf "%.4f%s\n", t, line }' conf$i
  i=$((i + 1))
done | awk '{ print; for (c = 2; c <= NF; c++) if (NR == 1 || $c < m[c]) m[c] = $c; n = NF }
  END { printf "# smallest:"; for (c = 2; c <= n; c++) printf " %.4f", m[c]; printf "\n" }'
