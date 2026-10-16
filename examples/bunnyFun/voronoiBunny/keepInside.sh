#!/bin/sh
# Keeps the cells whose centre is inside a shape (a cookie cutter), from the
# files written by tess2shape: the particles of a conf-file and their shapes.
# The cutter is read in a shape file, placed at the origin and scaled by a
# homothety; its faces must form a closed surface. A centre is inside when rays
# cast from it cross the faces an odd number of times (three rays vote, as in
# generatePacking:grid_inside).
# usage: sh keepInside.sh cutter.shp CutterName homothety cells.txt cells.shp out.txt out.shp
cutFile=$1; cutName=$2; h=$3; inCells=$4; inShapes=$5; outCells=$6; outShapes=$7
export LC_ALL=C

# 1. the names of the cells kept, and their lines
awk -v name="$cutName" -v h="$h" -v outCells="$outCells" '
function cross3(ax, ay, az, bx, by, bz) { cx = ay * bz - az * by; cy = az * bx - ax * bz; cz = ax * by - ay * bx }
function crossings(px, py, pz, dx, dy, dz,    f, k, a, b, c, e1x, e1y, e1z, e2x, e2y, e2z, det, inv, sx, sy, sz, u, v, qx, qy, qz, hx, hy, hz, t, n) {
  n = 0
  for (f = 0; f < nf; f++) {
    a = fa[f]; b = fb[f]; c = fc[f]
    e1x = X[b] - X[a]; e1y = Y[b] - Y[a]; e1z = Z[b] - Z[a]
    e2x = X[c] - X[a]; e2y = Y[c] - Y[a]; e2z = Z[c] - Z[a]
    cross3(dx, dy, dz, e2x, e2y, e2z); hx = cx; hy = cy; hz = cz
    det = e1x * hx + e1y * hy + e1z * hz
    if (det < 1e-18 && det > -1e-18) continue
    inv = 1 / det
    sx = px - X[a]; sy = py - Y[a]; sz = pz - Z[a]
    u = (sx * hx + sy * hy + sz * hz) * inv
    if (u < 0 || u > 1) continue
    cross3(sx, sy, sz, e1x, e1y, e1z); qx = cx; qy = cy; qz = cz
    v = (dx * qx + dy * qy + dz * qz) * inv
    if (v < 0 || u + v > 1) continue
    t = (e2x * qx + e2y * qy + e2z * qz) * inv
    if (t > 0) n++
  }
  return n
}
FNR == NR {                                   # first file: the cutter
  if ($1 == "name") inCut = ($2 == name)
  if (inCut && $1 == "nv") { nv = $2; for (i = 0; i < nv; i++) { getline; X[i] = $1 * h; Y[i] = $2 * h; Z[i] = $3 * h } }
  if (inCut && $1 == "nf") { m = $2; for (i = 0; i < m; i++) { getline; for (k = 3; k < $1 + 1; k++) { fa[nf] = $2; fb[nf] = $k; fc[nf] = $(k + 1); nf++ } } }
  next
}
/^Particles/ { n = $2; reading = 1; next }
reading && n > 0 {
  n--
  px = $5; py = $6; pz = $7; votes = 0
  if (crossings(px, py, pz, 0.8017, 0.4523, 0.3906) % 2 == 1) votes++
  if (crossings(px, py, pz, -0.3172, 0.8761, -0.3631) % 2 == 1) votes++
  if (crossings(px, py, pz, -0.5517, -0.3291, 0.7663) % 2 == 1) votes++
  if (votes >= 2) { kept[++nk] = $0; print $1 }
  total++
  next
}
END {
  print "Particles " nk > outCells
  for (i = 1; i <= nk; i++) print kept[i] > outCells
  printf "%d cells kept out of %d\n", nk, total > "/dev/stderr"
}' "$cutFile" "$inCells" > keep.names

# 2. their shapes
awk 'FNR == NR { keep[$1] = 1; next }
  /^</ { block = $0 "\n"; inside = 1; next }
  inside { block = block $0 "\n"; if ($1 == "name") name = $2 }
  /^>/ { if (name in keep) printf "%s", block; inside = 0 }' keep.names "$inShapes" > "$outShapes"
rm -f keep.names
