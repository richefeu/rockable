#!/bin/sh
# sheet.sh <work folder> <output png>
#
# Puts four frames (first, 1/3, 2/3, last) side by side, to check an
# animation at a glance before looking at the GIF.
O=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
cd "$1"
n=$(ls frame*.tga | wc -l); a=0; b=$((n / 3)); c=$((2 * n / 3)); d=$((n - 1))
magick $(printf "frame%04d.tga frame%04d.tga frame%04d.tga frame%04d.tga" $a $b $c $d) \
  -crop 768x708+0+0 +repage -resize 40% +append "$O"
