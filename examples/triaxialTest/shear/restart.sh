#!/bin/sh
# Builds input.txt, the initial state of the shear, from a conf-file of the
# compression (conf30 by default):
#
#   sh restart.sh [conf-file]
#
# Any conf-file is a valid input file of Rockable: it holds the parameters, the
# particles and their contacts. Only a few lines are changed:
#   - the time and the number of the conf-files start again from 0,
#   - the duration and the period of the conf-files of the shear,
#   - the friction between the spheres is switched on.
# drivingSystem.txt, shapes.txt and dataExtractors.txt are not in the
# conf-file: they are read in the folder of the computation, here.

CONF=${1:-../compression/conf30}

sed -e 's/^t .*/t 0/' \
    -e 's/^iconf .*/iconf 0/' \
    -e 's/^tmax .*/tmax 1.0/' \
    -e 's/^interConf .*/interConf 0.02/' \
    -e 's/^muContact 1 1 .*/muContact 1 1 0.5/' \
    "$CONF" > input.txt

echo "input.txt written from $CONF"
