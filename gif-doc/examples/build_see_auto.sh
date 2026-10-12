#!/bin/sh
# build_see_auto.sh [build folder]
#
# Builds see_auto, a copy of see that renders the conf-files to TGA frames
# and exits (see autoshot.inc). see.cpp is not modified: a copy is made in
# this folder, and compiled with the flags of the target `see` of the build
# folder (default: ../../BUILD, made by install_rockable.sh).
set -e
H=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$H/../.." && pwd)
B=$(cd "${1:-$R/BUILD}" && pwd)
D=$B/CMakeFiles/see.dir
[ -f "$D/flags.make" ] || { echo "no $D/flags.make: build Rockable with see first"; exit 1; }
[ -f "$B/libRockable_core.a" ] || { echo "no $B/libRockable_core.a"; exit 1; }

# autoShot is inserted before main, and armed just before glutMainLoop
awk -v inc="$H/autoshot.inc" '
  /^int main\(/ { while ((getline l < inc) > 0) print l }
  /^ *glutMainLoop\(\);/ { print "  glutTimerFunc(800, autoShot, 0);" }
  { print }' "$R/src/Apps/see.cpp" > "$H/see_auto.cpp"
grep -q "glutTimerFunc(800, autoShot" "$H/see_auto.cpp" || { echo "insertion failed"; exit 1; }

CXX=$(sed -n 's/^# compile CXX with //p' "$D/flags.make")
get() { sed -n "s/^$1 = //p" "$D/flags.make"; }
DEFS=$(get CXX_DEFINES); INCS=$(get CXX_INCLUDES); FLAGS=$(get CXX_FLAGS)
eval "$CXX $DEFS -I$R/src/Apps $INCS $FLAGS -c \"$H/see_auto.cpp\" -o \"$H/see_auto.o\""

# same link line as see, with see_auto.o and libRockable_core.a of the build folder
LINK=$(sed -e "s#CMakeFiles/see.dir/src/Apps/see.cpp.o#$H/see_auto.o#" -e "s#-o see #-o $H/see_auto #" \
           -e "s# libRockable_core.a# $B/libRockable_core.a#" "$D/link.txt")
eval "$LINK"
rm -f "$H/see_auto.o"
echo "built $H/see_auto"
