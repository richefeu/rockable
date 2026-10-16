#!/bin/sh
# Profile of the deposit in a conf-file: the height of the free surface along
# the channel, by slices of 1 cm in x.
#
#   sh profile.sh [conf-file]          (conf65 by default)
#
# The height of a slice is the highest centre of a grain in the slice, plus the
# radius of a grain (about 5 mm). The front of the deposit is the end of the
# first empty slice: the few grains that rolled further away are left out.

# Numbers are read and written with a decimal point, whatever the language settings
export LC_ALL=C

CONF=${1:-conf65}

awk '
  $1 == "Grain" {
    k = int($5 / 0.01)
    if (!(k in h) || $6 > h[k]) h[k] = $6
  }
  END {
    print "# x (m)  height (m)"
    for (k = 0; k in h; k++) printf "%.3f %.4f\n", (k + 0.5) * 0.01, h[k] + 0.005
    printf "# front of the deposit: x = %.3f m\n", k * 0.01
  }' "$CONF"
