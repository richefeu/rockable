#!/bin/sh
# render.sh <work folder> <see.json> <output gif>
#
# Renders the conf-files of the work folder with see_auto (a window opens
# for a while), then makes a light GIF with ffmpeg and gifsicle.
# Environment:
#   STEP, LAST     conf-files between two frames (1), last conf-file (all)
#   CMIN, CMAX     fixed colour range (colorMode of see.json)
#   WIDTH, FPS     width of the GIF in pixels (300), frames per second (12)
#   FRAMES         number of frames kept (all)
#   COLORS, LOSSY  palette size (32), lossy level of gifsicle (60)
#   DITHER         ffmpeg dithering (bayer:bayer_scale=4, or none)
#   FILTER         extra ffmpeg filters placed before the scaling
set -e
H=$(cd "$(dirname "$0")" && pwd)
J=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
O=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
[ -x "$H/see_auto" ] || { echo "run build_see_auto.sh first"; exit 1; }
cd "$1"
cp "$J" see.json
rm -f frame*.tga
export AUTOSHOT_STEP=${STEP:-1}
[ -n "$LAST" ] && export AUTOSHOT_LAST=$LAST
[ -n "$CMIN" ] && [ -n "$CMAX" ] && export AUTOSHOT_CMIN=$CMIN AUTOSHOT_CMAX=$CMAX
"$H/see_auto" conf0 > see.log 2>&1
# the frames are 768x768: the bottom strip (text zone) is dropped
F="crop=iw:ih-60:0:0${FILTER:+,$FILTER},scale=${WIDTH:-300}:-1:flags=lanczos"
ffmpeg -loglevel error -y -framerate "${FPS:-12}" -i frame%04d.tga ${FRAMES:+-frames:v $FRAMES} \
  -vf "$F,split[a][b];[a]palettegen=max_colors=${COLORS:-32}:stats_mode=diff[p];[b][p]paletteuse=dither=${DITHER:-bayer:bayer_scale=4}" tmp.gif
gifsicle -O3 --lossy="${LOSSY:-60}" --loopcount=forever tmp.gif -o "$O"
rm -f tmp.gif
echo "  $(basename "$O"): $(ls frame*.tga | wc -l | tr -d ' ') frames, $(($(wc -c < "$O") / 1024)) KB"
