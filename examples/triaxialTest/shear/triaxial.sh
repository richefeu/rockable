#!/bin/sh
# The results of the triaxial test, from SampleState.txt (written by
# 'postpro state.txt'):
#
#   sh triaxial.sh > triaxial.txt
#
# Columns of the output:
#   1 time
#   2 axial strain        eps_y = 1 - Ly / Ly0   (positive in compression)
#   3 volumetric strain   eps_v = 1 - V / V0     (negative when the sample dilates)
#   4 mean stress p, in Pa
#   5 deviatoric stress  syy - (sxx + szz) / 2, in Pa
#   6 stress ratio q / p
#   7 solid fraction
#   8 coordination number Z

# Numbers are read and written with a decimal point, whatever the language settings
export LC_ALL=C

# In SampleState.txt: 3 Ly, 5 V, 6 solid fraction, 7 Z, 9 p, 10 q,
# 11 sxx, 12 syy, 13 szz (the lines starting with # are comments)
awk '
  BEGIN { print "# time  eps_y  eps_v  p  syy-(sxx+szz)/2  q/p  solid_fraction  Z" }
  /^#/ { next }
  Ly0 == "" { Ly0 = $3; V0 = $5 }
  {
    printf "%.3f %.4f %.4f %.1f %.1f %.4f %.4f %.3f\n",
           $1, 1 - $3 / Ly0, 1 - $5 / V0, $9, $12 - 0.5 * ($11 + $13), $10 / $9, $6, $7
  }' SampleState.txt
