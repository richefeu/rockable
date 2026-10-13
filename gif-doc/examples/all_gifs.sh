#!/bin/sh
# all_gifs.sh [name ...]
#
# Makes the animation.gif of the examples (all of them, or only the ones
# named, e.g. `sh all_gifs.sh helloworld driving`). Each example is run in
# work/<name> (a shortened run when the full one is long), rendered with
# see_auto, and the GIF is written in work/<name>.gif together with a contact
# sheet work/<name>.png. Copy the GIFs into the examples once checked:
#   sh all_gifs.sh install
#
# Needs: Rockable built in ../../BUILD, or in the folder given by the variable BIN
# (with see), see_auto (build_see_auto.sh),
# ffmpeg, gifsicle, ImageMagick (magick).
# Times on a laptop: a few seconds to a few minutes per example, except
# softBallImpactVoronoi (~4 min with 4 threads), DominoFun (~12 min),
# prepro_cpptools (~10 min) and PerioPoly (~15 min).
set -e
H=$(cd "$(dirname "$0")" && pwd)
X=$(cd "$H/../../examples" && pwd)
W=$H/work
BIN=$(cd "${BIN:-$H/../../BUILD}" && pwd)
export ROCKABLE=${ROCKABLE:-$BIN/rockable}
mkdir -p "$W"

# camera: json <name> '<center>' '<eye>' [extra json]
json() {
  printf '{"window": {"width": 768, "height": 768}, "camera": {"center": [%s], "eye": [%s], "up": [0, 1, 0]}, "show_background": 0%s}\n' \
    "$2" "$3" "${4:+, $4}" > "$W/$1.json"
}
# standard case: run, then render, then contact sheet
# std <name> <example folder> <interConf> <tmax or ''>   (env of run_example.sh and render.sh)
std() {
  sh "$H/run_example.sh" "$2" "$W/$1" "$3" "$4"
  sh "$H/render.sh" "$W/$1" "$W/$1.json" "$W/$1.gif"
  sh "$H/sheet.sh" "$W/$1" "$W/$1.png"
}

helloworld() {
  json helloworld '-0.15, 0.22, 0' '0.35, 0.75, 2.4'
  WIDTH=360 FPS=15 COLORS=64 LOSSY=40 std helloworld "$X/helloworld" 0.01 ''
}
helloworldPlus() {
  json helloworldPlus '0.25, 0.3, 0' '0.75, 0.9, 3.0'
  WIDTH=320 FPS=15 COLORS=32 LOSSY=80 std helloworldPlus "$X/helloworldPlus" 0.01 ''
}
driving() {
  json driving '0, 0, 0' '0, 0.3, 6.0' '"alpha_fixparticles": 0.6'
  WIDTH=320 FPS=15 COLORS=32 LOSSY=60 std driving "$X/driving" 0.01 ''
}
isoComprSpheres() {
  json isoComprSpheres '0.025, 0.025, 0.025' '0.084, 0.0742, 0.1234' '"show_driven": 0'
  std isoComprSpheres "$X/isoComprSpheres" 0.005 0.25
}
prepro_cpptools() {
  # the input file and the shapes are written by generator (genAll.txt)
  json prepro_cpptools '0.125, 0.125, 0.125' '0.3541, 0.3312, 0.5832' '"show_driven": 0'
  PRE="$BIN/generator genAll.txt > gen.log 2>&1" JOBS=4 \
    sh "$H/run_example.sh" "$X/prepro_cpptools" "$W/prepro_cpptools" 0.01 0.3
  LAST=30 FRAMES=30 WIDTH=220 FPS=10 COLORS=12 DITHER=none LOSSY=150 \
    sh "$H/render.sh" "$W/prepro_cpptools" "$W/prepro_cpptools.json" "$W/prepro_cpptools.gif"
  sh "$H/sheet.sh" "$W/prepro_cpptools" "$W/prepro_cpptools.png"
}
DominoFun() {
  json DominoFun '0.36, 0.06, 0' '0.65, 0.62, 0.95' '"alpha_fixparticles": 0.5'
  # the input file and the shapes are written by generator (generate.txt)
  PRE="$BIN/generator generate.txt > gen.log 2>&1" \
    JOBS=4 COLORS=32 LOSSY=80 std DominoFun "$X/DominoFun" 0.025 2.5
}
stickVerticesInClusters() {
  json stickVerticesInClusters '0, 0.0085, 0' '0, 0.0181, 0.0640' '"alpha_fixparticles": 0.4'
  std stickVerticesInClusters "$X/stickVerticesInClusters" 0.0006 0.06
}
isoComprSpheresPeriodic() {
  json isoComprSpheresPeriodic '0.025, 0.025, 0.025' '0.0817, 0.0723, 0.1196' '"show_periodicBox": 1'
  std isoComprSpheresPeriodic "$X/isoComprSpheresPeriodic" 0.004 0.3
}
PerioPoly() {
  # the compression of the cell is over at t = 0.18 s
  json PerioPoly '3.7, 3.7, 3.7' '12.5, 11.5, 20.5' '"show_periodicBox": 1'
  JOBS=2 LAST=36 WIDTH=260 FPS=10 COLORS=24 DITHER=none LOSSY=100 std PerioPoly "$X/PerioPoly" 0.006 0.6
}
softBallImpact() {
  json softBallImpact '1, 1, 1' '2.35, 2.1, 5.4' '"colorMode": 4, "alpha_fixparticles": 0.5'
  LAST=45 CMIN=-0.4 CMAX=0.4 COLORS=32 LOSSY=80 std softBallImpact "$X/softParticles/softBallImpact" 0.01 ''
}
softBallImpactVoronoi() {
  json softBallImpactVoronoi '1, 1, 1' '2.35, 2.3, 5.6' '"colorMode": 4, "alpha_fixparticles": 0.5'
  JOBS=4 CMIN=-0.4 CMAX=0.4 COLORS=32 LOSSY=80 std softBallImpactVoronoi "$X/softParticles/softBallImpactVoronoi" 0.008 0.36
}
elasticCube() {
  json elasticCube '0, -0.5, 0' '2.6, -0.3, 3.4' '"colorMode": 5, "alpha_fixparticles": 0.25'
  CMIN=-4000 CMAX=0 WIDTH=280 COLORS=48 std elasticCube "$X/softParticles/elasticCube" 0.04 ''
}
viscousCube() {
  json viscousCube '0, -0.5, 0' '2.6, -0.3, 3.4' '"colorMode": 5, "alpha_fixparticles": 0.25'
  CMIN=-4000 CMAX=0 WIDTH=280 COLORS=48 std viscousCube "$X/softParticles/viscousCube" 0.04 ''
}
bondedBeam() {
  json bondedBeam '0.1, 0.0035, 0' '0.1590, 0.0773, 0.2952' '"colorMode": 5, "alpha_fixparticles": 0.5'
  CMIN=-5000 CMAX=5000 WIDTH=320 COLORS=48 std bondedBeam "$X/softParticles/bondedBeam" 0.006 ''
}
sweepable() {
  # four cases of the sweep (en2, mu), side by side
  json sweepable '0.25, 0.3, 0' '0.75, 0.9, 3.0'
  F=/System/Library/Fonts/Supplemental/Arial.ttf
  [ -f "$F" ] || F=$(fc-match -f '%{file}' sans 2>/dev/null || true)
  D=$W/sweepable; rm -rf "$D"; mkdir -p "$D/tiles"
  for c in "0.01 0.1" "0.01 0.5" "0.6 0.1" "0.6 0.5"; do
    set -- $c
    SED="s/^en2Contact 0 1 .*/en2Contact 0 1 $1/;s/^muContact 0 1 .*/muContact 0 1 $2/" \
      sh "$H/run_example.sh" "$X/sweepable" "$D/en2_$1_mu_$2" 0.02 ''
    (cd "$D/en2_$1_mu_$2" && cp "$W/sweepable.json" see.json && "$H/see_auto" conf0 > see.log 2>&1)
  done
  n=$(ls "$D/en2_0.01_mu_0.1" | grep -c '^frame')
  i=0
  while [ $i -lt "$n" ]; do
    f=$(printf frame%04d.tga $i)
    for d in en2_0.01_mu_0.1 en2_0.01_mu_0.5 en2_0.6_mu_0.1 en2_0.6_mu_0.5; do
      l=$(echo $d | sed 's/_mu_/   mu = /;s/en2_/en2 = /')
      magick "$D/$d/$f" -crop 768x708+0+0 +repage -resize 200x ${F:+-font "$F"} -gravity north -pointsize 14 \
        -annotate +0+6 "$l" -bordercolor gray80 -border 1 "$D/tiles/$d.png"
    done
    magick \( "$D/tiles/en2_0.01_mu_0.1.png" "$D/tiles/en2_0.01_mu_0.5.png" +append \) \
           \( "$D/tiles/en2_0.6_mu_0.1.png" "$D/tiles/en2_0.6_mu_0.5.png" +append \) -append "$D/tiles/$(printf m%04d.png $i)"
    i=$((i + 1))
  done
  ffmpeg -loglevel error -y -framerate 12 -i "$D/tiles/m%04d.png" \
    -vf "split[a][b];[a]palettegen=max_colors=32:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4" "$D/tmp.gif"
  gifsicle -O3 --lossy=80 --loopcount=forever "$D/tmp.gif" -o "$W/sweepable.gif"
  cp "$D/tiles/m0025.png" "$W/sweepable.png"
  echo "  sweepable.gif: $(($(wc -c < "$W/sweepable.gif") / 1024)) KB"
}
rockfall() {
  json rockfall '12, -1, 0' '12, 26, 72' '"alpha_fixparticles": 0.6'
  FILTER='crop=iw:400:0:150' WIDTH=560 COLORS=32 LOSSY=80 std rockfall "$X/rockfall" 0.1 10
}
fragmentation() {
  json fragmentation '0, 0.03, 0' '0.16, 0.13, 0.29' '"alpha_fixparticles": 0.5'
  WIDTH=260 COLORS=24 LOSSY=100 std fragmentation "$X/fragmentation" 0.00025 ''
}
columnCollapse() {
  # from the settled column (conf24), without the walls, in a band
  json columnCollapse '0.14, 0.05, 0.025' '0.14, 0.11, 0.42' '"show_driven": 0'
  FIRST=24 FILTER='crop=iw:320:0:190' WIDTH=440 COLORS=24 DITHER=none LOSSY=120 std columnCollapse "$X/columnCollapse" '' ''
}
bunnyPile() {
  # the bunnies fall into the tube, settle, and spread when the tube is lifted (the tube and the floor
  # are hidden: see draws them in the depth buffer, and they would hide the bunnies inside); with skins
  json bunnyPile '0, 0.1, 0' '0.3, 0.3, 0.55' '"show_driven": 0'
  PRE="$BIN/stl2shape -r 0.001 -c $X/bunnyFun/bunny.stl > stl2shape.log && sed 's/^name .*/name Bunny/' mesh.shp > bunny.shp && $BIN/shapeSurvey bunny.shp -P bunny.shp > ss.log && $BIN/generator genShapes.txt > gen.log 2>&1 && cat walls.shp bunny.shp > shapes.shp && $BIN/generator genPile.txt >> gen.log 2>&1 && $BIN/shape2mesh shapes.shp --rmsh -s > s2m.log" \
    JOBS=4 sh "$H/run_example.sh" "$X/bunnyFun/bunnyPile" "$W/bunnyPile" '' ''
  FIRST=10 WIDTH=300 COLORS=32 LOSSY=80 sh "$H/render.sh" "$W/bunnyPile" "$W/bunnyPile.json" "$W/bunnyPile.gif"
  sh "$H/sheet.sh" "$W/bunnyPile" "$W/bunnyPile.png"
}
postProcessing() {
  # the force chains: tubes between the centres, as thick as the normal force, in translucent spheres
  json postProcessing '0.017, 0.017, 0.017' '0.053, 0.047, 0.077' '"show_driven": 0, "show_normal_forces": 1, "alpha_particles": 0.15'
  FIRST=16 STEP=2 WIDTH=260 COLORS=32 LOSSY=80 std postProcessing "$X/postProcessing" '' ''
}
triaxialTest() {
  # the shear only: the compression is run first, and restart.sh builds the input file of the shear
  json triaxialTest '0.026, 0.022, 0.026' '0.094, 0.056, 0.126' '"show_driven": 0'
  sh "$H/run_example.sh" "$X/triaxialTest/compression" "$W/triaxialTest_compression" '' ''
  PRE="sh restart.sh $W/triaxialTest_compression/conf30" STEP=3 FPS=8 WIDTH=240 COLORS=24 DITHER=none LOSSY=120 std triaxialTest "$X/triaxialTest/shear" '' ''
}
complexShape() {
  # the rock of rock.stl spinning freely (the input.txt of the example), drawn twice side by side:
  # with spheres, tubes and slabs (left), and with its skin mesh, rock.rmsh (right)
  json complexShape '0, 0, 0' '0, 0.03, 0.1'
  D=$W/complexShape
  PRE="$BIN/stl2shape -r 0.5 -c rock.stl > stl2shape.log && mv mesh.shp rock.shp && $BIN/shapeSurvey rock.shp -P rock.shp > ss.log" \
    sh "$H/run_example.sh" "$X/complexShape" "$D" '' ''
  F=/System/Library/Fonts/Supplemental/Arial.ttf
  [ -f "$F" ] || F=$(fc-match -f '%{file}' sans 2>/dev/null || true)
  for v in primitives skin; do
    rm -rf "$D/$v"; mkdir -p "$D/$v"
    (cd "$D" && cp conf* rock.shp "$v/")
    [ $v = skin ] && (cd "$D/$v" && "$BIN/shape2mesh" rock.shp --rmsh > shape2mesh.log)
    (cd "$D/$v" && cp "$W/complexShape.json" see.json && AUTOSHOT_STEP=4 AUTOSHOT_LAST=49 "$H/see_auto" conf0 > see.log 2>&1)
  done
  mkdir -p "$D/tiles"; rm -f "$D"/tiles/*
  n=$(ls "$D/skin" | grep -c '^frame')
  i=0
  while [ $i -lt "$n" ]; do
    f=$(printf frame%04d.tga $i)
    for v in primitives skin; do
      l="spheres, tubes, slabs"; [ $v = skin ] && l="skin (rock.rmsh)"
      magick "$D/$v/$f" -crop 768x708+0+0 +repage -gravity center -crop 52%x56%+0+0 +repage -resize 230x \
        ${F:+-font "$F"} -gravity north -pointsize 14 -fill black -undercolor "#ffffffb0" -annotate +0+4 " $l " "$D/tiles/$v.png"
    done
    magick "$D/tiles/primitives.png" "$D/tiles/skin.png" +append "$D/tiles/$(printf m%04d.png $i)"
    i=$((i + 1))
  done
  ffmpeg -loglevel error -y -framerate 4 -i "$D/tiles/m%04d.png" \
    -vf "split[a][b];[a]palettegen=max_colors=24:stats_mode=diff[p];[b][p]paletteuse=dither=none" "$D/tmp.gif"
  gifsicle -O3 --lossy=200 --loopcount=forever "$D/tmp.gif" -o "$W/complexShape.gif"
  cp "$D/tiles/m0006.png" "$W/complexShape.png"
  echo "  complexShape.gif: $n frames, $(($(wc -c < "$W/complexShape.gif") / 1024)) KB"
}

# where each GIF goes in the examples
dest() {
  case $1 in
    softBallImpact | softBallImpactVoronoi | elasticCube | viscousCube | bondedBeam) echo "$X/softParticles/$1" ;;
    bunnyPile) echo "$X/bunnyFun/$1" ;;
    *) echo "$X/$1" ;;
  esac
}
ALL="helloworld helloworldPlus driving isoComprSpheres prepro_cpptools DominoFun stickVerticesInClusters
     isoComprSpheresPeriodic PerioPoly softBallImpact softBallImpactVoronoi bondedBeam elasticCube viscousCube
     sweepable complexShape rockfall triaxialTest postProcessing fragmentation columnCollapse bunnyPile"

if [ "$1" = install ]; then
  for e in $ALL; do
    [ -f "$W/$e.gif" ] && cp "$W/$e.gif" "$(dest $e)/animation.gif" && echo "installed $e"
  done
  exit 0
fi
for e in ${*:-$ALL}; do
  echo "== $e"
  ( $e )
done
