#!/bin/sh
# Distribution of the normal contact forces in a conf-file, normalised by their
# mean: the probability density of fn/<fn>, by classes of width 0.25.
#
#   sh forceDistribution.sh [conf-file]     (conf60 by default)
#
# It reads the block 'Interactions' of the conf-file: each line is one
# interaction, and the normal force fn is its 16th column (see the page
# syntaxConf of the documentation). Only the contacts (dn < 0, 9th column)
# with a repulsive force are kept.

export LC_ALL=C
CONF=${1:-conf60}

awk -v w=0.25 '
  $1 == "Interactions" { n = $2; for (k = 0; k < n; k++) { getline; if ($9 < 0 && $16 > 0) fn[N++] = $16 } }
  END {
    if (N == 0) { print "# no contact"; exit }
    for (k = 0; k < N; k++) sum += fn[k]
    mean = sum / N
    for (k = 0; k < N; k++) { c = int(fn[k] / mean / w); count[c]++; if (c > cmax) cmax = c }
    printf "# %d contacts, mean normal force %g\n", N, mean
    print "# fn/<fn> (centre of the class)  probability density"
    for (c = 0; c <= cmax; c++) printf "%.3f %.4f\n", (c + 0.5) * w, count[c] / N / w
  }' "$CONF"
