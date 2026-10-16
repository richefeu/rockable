#!/bin/sh
# Makes a hollow shape with a hole from a shape written by stl2shape:
#  - the faces whose centre is closer than R to the point (X, Y, Z) are removed,
#    and so are the edges and the vertices that no longer belong to any face;
#  - the shape is marked isSurface: it is a shell, of thickness twice its radius;
#  - it is marked preCompDone y, with an OBB that bounds it, the volume of the
#    shell (the area of its faces times its thickness), and the inertia of a
#    box (the shape is meant for a driven body, whose mass hardly matters).
# usage: sh hollow.sh in.shp out.shp X Y Z R
in=$1; out=$2
export LC_ALL=C
awk -v X="$3" -v Y="$4" -v Z="$5" -v R="$6" '
$1 == "name"   { name = $2 }
$1 == "radius" { radius = $2 }
$1 == "nv" { nv = $2; for (i = 0; i < nv; i++) { getline; x[i] = $1; y[i] = $2; z[i] = $3 } }
$1 == "ne" { ne = $2; for (i = 0; i < ne; i++) { getline; ea[i] = $1; eb[i] = $2 } }
$1 == "nf" {
  nf = $2
  for (f = 0; f < nf; f++) {
    getline
    cx = 0; cy = 0; cz = 0
    for (k = 2; k <= NF; k++) { cx += x[$k]; cy += y[$k]; cz += z[$k] }
    k = NF - 1; cx /= k; cy /= k; cz /= k
    if ((cx - X)^2 + (cy - Y)^2 + (cz - Z)^2 < R * R) { removed++; continue }
    for (k = 3; k < NF; k++) {                 # area of the face, as a fan of triangles
      ux = x[$k] - x[$2]; uy = y[$k] - y[$2]; uz = z[$k] - z[$2]
      vx = x[$(k + 1)] - x[$2]; vy = y[$(k + 1)] - y[$2]; vz = z[$(k + 1)] - z[$2]
      area += 0.5 * sqrt((uy * vz - uz * vy)^2 + (uz * vx - ux * vz)^2 + (ux * vy - uy * vx)^2)
    }
    face[kept++] = $0
    for (k = 2; k <= NF; k++) {
      used[$k] = 1
      a = $k; b = (k < NF) ? $(k + 1) : $2      # the edges of the face
      onface[(a < b) ? a " " b : b " " a] = 1
    }
  }
}
END {
  for (i = 0; i < nv; i++) if (i in used) { newid[i] = n++; if (n == 1) { x0 = x1 = x[i]; y0 = y1 = y[i]; z0 = z1 = z[i] }
    if (x[i] < x0) x0 = x[i]; if (x[i] > x1) x1 = x[i]; if (y[i] < y0) y0 = y[i]; if (y[i] > y1) y1 = y[i]
    if (z[i] < z0) z0 = z[i]; if (z[i] > z1) z1 = z[i] }
  for (i = 0; i < ne; i++) { a = ea[i]; b = eb[i]; if (((a < b) ? a " " b : b " " a) in onface) edge[m++] = newid[a] " " newid[b] }
  print "<"
  print "name " name
  print "radius " radius
  print "isSurface"
  print "preCompDone y"
  print "nv " n
  for (i = 0; i < nv; i++) if (i in used) print x[i], y[i], z[i]
  print "ne " m
  for (i = 0; i < m; i++) print edge[i]
  print "nf " kept
  for (f = 0; f < kept; f++) {
    split(face[f], w, " "); line = w[1]
    for (k = 2; k <= w[1] + 1; k++) line = line " " newid[w[k]]
    print line
  }
  printf "obb.extent %g %g %g\n", (x1 - x0) / 2 + radius, (y1 - y0) / 2 + radius, (z1 - z0) / 2 + radius
  print "obb.e1 1 0 0"
  print "obb.e2 0 1 0"
  print "obb.e3 0 0 1"
  printf "obb.center %g %g %g\n", (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2
  printf "volume %g\n", area * 2 * radius
  a = x1 - x0; b = y1 - y0; c = z1 - z0
  printf "I/m %g %g %g\n", (b * b + c * c) / 12, (a * a + c * c) / 12, (a * a + b * b) / 12
  print ">"
  printf "%d faces removed, %d kept; %d vertices, %d edges\n", removed, kept, n, m > "/dev/stderr"
}' "$in" > "$out"
