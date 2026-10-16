#!/bin/sh
# Builds input.txt and shapes.txt: the parameters (params.txt), the plate, the
# ball, and the cells kept in the bunny (kept.txt, kept.shp), lifted by 5.15 cm
# so that the lowest ones lie on the plate. The cells, all in the cluster 0, are
# glued by stickBCM.
# usage: sh makeInput.sh
export LC_ALL=C
cat others.shp kept.shp > shapes.txt
n=$(awk '/^Particles/ { print $2 }' kept.txt)
{
  cat params.txt
  echo
  echo "Particles $((n + 2))"
  echo "Platen 0 9 1  0 -0.002 0  0 0 0  0 0 0  1 0 0 0  0 0 0  0 0 0"
  echo "Ball 2 1 1  -0.075 0.05 0.2  0 0 -60  0 0 0  1 0 0 0  0 0 0  0 0 0"
  awk 'NR > 1 { $6 = sprintf("%.10g", $6 + 0.0515); print }' kept.txt
  echo
  echo "# The cells of the same cluster are glued, face to face"
  echo "stickBCM 1e-4"
} > input.txt
echo "input.txt written: $n cells, the plate and the ball"
