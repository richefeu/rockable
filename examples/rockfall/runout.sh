#!/bin/sh
# Run-out of the blocks: for each file written by TrackBody, the largest speed,
# the time when the block stops (speed below 0.05 m/s for good), and the
# distance it travelled beyond the foot of the slope (x = 15.32 m).
#
#   sh runout.sh [files]      (sphere.txt cube.txt slab.txt by default)
#
# Columns of a TrackBody file: 1 t, 2-4 position, 5-7 velocity, ...

export LC_ALL=C
FILES=${*:-"sphere.txt cube.txt slab.txt"}
echo "# block  max speed (m/s)  stop time (s)  run-out beyond the foot (m)"
for f in $FILES; do
  awk -v name="${f%.txt}" -v xfoot=15.32 '
    { v = sqrt($5 * $5 + $6 * $6 + $7 * $7); if (v > vmax) vmax = v
      if (v > 0.05) tstop = -1; else if (tstop < 0) tstop = $1
      x = $2 }
    END { printf "%-8s %8.2f %14s %18.2f\n", name, vmax, (tstop < 0 ? "moving" : sprintf("%.1f", tstop)), x - xfoot }' "$f"
done
