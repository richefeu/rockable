#!/bin/sh
# The grains (the particles of group 0) in a conf-file: their number, the
# highest and the mean height of their centres, their mean distance to the
# axis of the tube, and the solid fraction of the tube up to the highest centre
# (meaningful while the tube is in place). The volume of each shape is read
# in shapes.shp.
# usage: sh pile.sh conf35
conf=${1:-conf35}
export LC_ALL=C
awk -v R=0.06 '
FNR == NR {                                  # first file: the shapes
  if ($1 == "name") name = $2
  if ($1 == "volume") V[name] = $2
  next
}
/^t / { t = $2 }
/^Particles/ { n = $2; reading = 1; next }
reading && /^#/ { next }
reading && n > 0 {
  n--
  if ($2 == 0) {                             # group 0: the grains
    ng++
    h = $4
    vol += V[$1] * h * h * h
    if ($6 > top) top = $6
    ymean += $6
    rmean += sqrt($5 * $5 + $7 * $7)
  }
  if (n == 0) reading = 0
}
END {
  printf "time                      %g s\n", t
  printf "grains                    %d\n", ng
  printf "highest centre            %.3f m\n", top
  printf "mean height of the centres %.3f m\n", ymean / ng
  printf "mean distance to the axis %.3f m\n", rmean / ng
  printf "solid fraction in the tube %.3f\n", vol / (3.14159265 * R * R * top)
}' shapes.shp "$conf"
