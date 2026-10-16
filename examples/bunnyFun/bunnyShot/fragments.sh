#!/bin/sh
# The fragments of the bunny in a conf-file: the cells (group 1) that are still
# linked by an interface belong to the same fragment. Gives the number of
# fragments, the number of single cells, and the sizes of the largest ones.
# usage: sh fragments.sh conf40
conf=${1:-conf40}
export LC_ALL=C
awk '
function find(i) { while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i] } return i }
/^Particles/ { n = $2; k = 0; reading = 1; next }
reading && /^#/ { next }
reading && n > 0 {
  if ($2 == 1) { cell[k] = 1; parent[k] = k; ncells++ }
  k++; n--
  if (n == 0) reading = 0
  next
}
/^Interfaces/ { m = $2; for (l = 0; l < m; l++) { getline; if (($1 in cell) && ($2 in cell)) { a = find($1); b = find($2); if (a != b) parent[a] = b } } }
END {
  for (i in cell) size[find(i)]++
  for (r in size) { nf++; if (size[r] == 1) single++; s[nf] = size[r] }
  # the sizes, largest first (a simple selection: there are few large ones)
  for (i = 1; i <= nf; i++) for (j = i + 1; j <= nf; j++) if (s[j] > s[i]) { t = s[i]; s[i] = s[j]; s[j] = t }
  printf "cells %d, fragments %d, single cells %d\n", ncells, nf, single
  printf "the largest fragments (cells):"
  for (i = 1; i <= nf && i <= 6; i++) printf " %d", s[i]
  printf "\n"
}' "$conf"
