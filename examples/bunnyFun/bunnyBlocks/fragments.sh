#!/bin/sh
# The fragments of the bunny in a conf-file: the blocks (group 1) that are still
# linked by an interface belong to the same fragment. The blocks do not have the
# same size: each fragment is given by its volume, in % of the bunny (the volumes
# of the shapes are read in shapes.txt).
# usage: sh fragments.sh [conf-file]          (conf60 by default)
conf=${1:-conf60}
export LC_ALL=C
awk '
function find(i) { while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i] } return i }
FNR == NR { if ($1 == "name") nm = $2; if ($1 == "volume") V[nm] = $2; next }
/^Particles/ { n = $2; k = 0; reading = 1; next }
reading && /^#/ { next }
reading && n > 0 {
  if ($2 == 1) { block[k] = 1; vol[k] = V[$1]; parent[k] = k; nblocks++; total += V[$1] }
  k++; n--
  if (n == 0) reading = 0
  next
}
/^Interfaces/ { m = $2; for (l = 0; l < m; l++) { getline; if (($1 in block) && ($2 in block)) { a = find($1); b = find($2); if (a != b) parent[a] = b } } }
END {
  for (i in block) { r = find(i); size[r] += vol[i]; nb[r]++ }
  for (r in size) { nf++; if (nb[r] == 1) single++; s[nf] = size[r]; c[nf] = nb[r] }
  # largest first (a simple selection: there are few fragments)
  for (i = 1; i <= nf; i++) for (j = i + 1; j <= nf; j++) if (s[j] > s[i]) { t = s[i]; s[i] = s[j]; s[j] = t; t = c[i]; c[i] = c[j]; c[j] = t }
  printf "blocks %d, fragments %d, single blocks %d\n", nblocks, nf, single
  printf "the largest fragments (%% of the volume, blocks):"
  for (i = 1; i <= nf && i <= 5; i++) printf "  %.1f (%d)", 100 * s[i] / total, c[i]
  printf "\n"
}' shapes.txt "$conf"
