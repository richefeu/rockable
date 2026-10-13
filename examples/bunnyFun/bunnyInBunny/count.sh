#!/bin/sh
# The small bunnies (group 0) inside and outside the big one (BigBunny), in a
# conf-file. A vertical ray is cast downwards from the centre of each small
# bunny: it crosses the shell of the big bunny an odd number of times when the
# small bunny is inside. The hole is on the back, above the small bunnies, so
# it does not change the count. A small bunny above the big one is outside.
# The positions are taken in the frame of the big bunny, which may move and turn.
# The angle printed is the angle of its rotation, between 0 and 180 degrees.
# usage: sh count.sh conf95
conf=${1:-conf95}
export LC_ALL=C
awk '
FNR == NR {                                   # first file: the shape BigBunny
  if ($1 == "name") inBig = ($2 == "BigBunny")
  if (inBig && $1 == "nv") { nv = $2; for (i = 0; i < nv; i++) { getline; x[i] = $1; y[i] = $2; z[i] = $3; if (i == 0 || y[i] > ymax) ymax = y[i] } }
  if (inBig && $1 == "nf") { nf = $2; for (f = 0; f < nf; f++) { getline; a[f] = $2; b[f] = $3; c[f] = $4 } }
  next
}
/^t / { t = $2 }
/^Particles/ { n = $2; reading = 1; next }
reading && /^#/ { next }
reading && n > 0 {
  n--
  if ($1 == "BigBunny") { h = $4; px = $5; py = $6; pz = $7; qw = $14; qx = -$15; qy = -$16; qz = -$17 }
  if ($2 == 0) { m++; sx[m] = $5; sy[m] = $6; sz[m] = $7 }
  if (n == 0) reading = 0
}
END {
  for (k = 1; k <= m; k++) {
    # the centre in the frame of the big bunny: v = conj(Q) (p - pos) / h
    vx = (sx[k] - px) / h; vy = (sy[k] - py) / h; vz = (sz[k] - pz) / h
    tx = 2 * (qy * vz - qz * vy); ty = 2 * (qz * vx - qx * vz); tz = 2 * (qx * vy - qy * vx)
    lx = vx + qw * tx + (qy * tz - qz * ty)
    ly = vy + qw * ty + (qz * tx - qx * tz)
    lz = vz + qw * tz + (qx * ty - qy * tx)
    cross = 0
    for (f = 0; f < nf; f++) {
      # does the line through (lx, lz), along y, cross the triangle below ly?
      x1 = x[a[f]]; z1 = z[a[f]]; x2 = x[b[f]]; z2 = z[b[f]]; x3 = x[c[f]]; z3 = z[c[f]]
      d = (z2 - z3) * (x1 - x3) + (x3 - x2) * (z1 - z3)
      if (d == 0) continue
      l1 = ((z2 - z3) * (lx - x3) + (x3 - x2) * (lz - z3)) / d
      l2 = ((z3 - z1) * (lx - x3) + (x1 - x3) * (lz - z3)) / d
      l3 = 1 - l1 - l2
      if (l1 < 0 || l2 < 0 || l3 < 0) continue
      if (l1 * y[a[f]] + l2 * y[b[f]] + l3 * y[c[f]] < ly) cross++
    }
    if (cross % 2 == 1 && ly < ymax) inside++
  }
  angle = 2 * atan2(sqrt(qx * qx + qy * qy + qz * qz), (qw < 0) ? -qw : qw) * 180 / 3.14159265
  printf "t = %.2f s, big bunny turned by %3.0f degrees: %d small bunnies inside, %d outside\n", t, angle, inside, m - inside
}' shapes.shp "$conf"
