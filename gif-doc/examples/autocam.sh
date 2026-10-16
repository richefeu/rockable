#!/bin/sh
# autocam.sh <conf-file> <dist> <dx> <dy> <dz> [extra json]  > see.json
#
# Writes a see.json whose camera looks at the centre of the bounding box of
# the particle centres, from the direction (dx, dy, dz), at a distance of
# dist times the size of the box. The extra json (for instance
# '"colorMode": 4, "show_driven": 0') is added to the parameters.
export LC_ALL=C
awk -v f="$2" -v dx="$3" -v dy="$4" -v dz="$5" -v extra="$6" '
  /^Particles/ { n = $2; reading = 1; next }
  reading && /^#/ { next }
  reading && k < n {
    for (i = 0; i < 3; i++) {
      x = $(5 + i)
      if (k == 0 || x < lo[i]) lo[i] = x
      if (k == 0 || x > hi[i]) hi[i] = x
    }
    k++
  }
  END {
    size = 0
    for (i = 0; i < 3; i++) { c[i] = (lo[i] + hi[i]) / 2; if (hi[i] - lo[i] > size) size = hi[i] - lo[i] }
    nd = sqrt(dx * dx + dy * dy + dz * dz)
    e[0] = c[0] + f * size * dx / nd; e[1] = c[1] + f * size * dy / nd; e[2] = c[2] + f * size * dz / nd
    printf "{\"window\": {\"width\": 768, \"height\": 768}, "
    printf "\"camera\": {\"center\": [%g, %g, %g], \"eye\": [%g, %g, %g], \"up\": [0, 1, 0]}, ", c[0], c[1], c[2], e[0], e[1], e[2]
    printf "\"show_background\": 0%s}\n", (extra != "" ? ", " extra : "")
  }' "$1"
