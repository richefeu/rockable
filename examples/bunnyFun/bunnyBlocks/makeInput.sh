#!/bin/sh
# Builds input.txt and shapes.txt: the parameters (params.txt), the plate, the
# cannonball, and the blocks cut by stl2blocks (blocks.txt, blocks.shp), glued
# by stickBCM. The cannonball (radius 6 cm) rolls without sliding towards -z, at
# V m/s: its angular velocity is -V/R about x. It starts 4 cm in front of the bunny.
# usage: sh makeInput.sh [V]   (V = 3 by default)
export LC_ALL=C
V=${1:-3}
cat others.shp blocks.shp > shapes.txt
n=$(awk '/^Particles/ { print $2 }' blocks.txt)
{
  cat params.txt
  echo
  echo "Particles $((n + 2))"
  echo "Platen 0 9 1  0 -0.002 0  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0"
  awk -v V=$V 'BEGIN { printf "Cannonball 2 1 1  -0.03 0.06 0.2  0 0 %g  0 0 0  1 0 0 0  %g 0 0  0 0 0\n", -V, -V / 0.06 }'
  awk 'NR > 3' blocks.txt
  echo
  echo "# The blocks of the same cluster are glued, face to face"
  echo "stickBCM 1e-4"
} > input.txt
echo "input.txt written: $n blocks, the plate and the cannonball at $V m/s"
