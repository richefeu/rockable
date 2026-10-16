#!/bin/sh
# run_example.sh <example folder> <work folder> [interConf] [tmax]
#
# Copies the example into the work folder (deleted first), removes its data
# extractors, sets interConf and tmax in the input file when they are given,
# and runs rockable there (output in run.log).
# Environment: INPUT (input file, default input.txt), PRE (a command run in
# the work folder first, for instance to generate the input file), SED (an
# extra sed script applied to the input file), JOBS (threads), ROCKABLE
# (executable, default ../../BUILD/rockable).
set -e
H=$(cd "$(dirname "$0")" && pwd)
ROCKABLE=${ROCKABLE:-$H/../../BUILD/rockable}
INPUT=${INPUT:-input.txt}
EX=$(cd "$1" && pwd); W=$2
rm -rf "$W"; cp -r "$EX" "$W"; cd "$W"
rm -f dataExtractors.txt conf* animation.gif
[ -n "$PRE" ] && eval "$PRE"
[ -n "$3" ] && sed -i.bak "s/^interConf .*/interConf $3/" "$INPUT"
[ -n "$4" ] && sed -i.bak "s/^tmax .*/tmax $4/" "$INPUT"
[ -n "$SED" ] && sed -i.bak "$SED" "$INPUT"
rm -f "$INPUT.bak"
echo "running $(basename "$EX") in $W"
"$ROCKABLE" "$INPUT" ${JOBS:+-j $JOBS} > run.log 2>&1
echo "  $(ls | grep -c '^conf[0-9]') conf-files"
