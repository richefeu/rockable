#!/usr/bin/env python3
"""Solid fraction of the sample of isoComprSpheres, during the compression.

It reads the positions of the three moving walls in xmax.txt, ymax.txt and
zmax.txt (written by the data extractors, see dataExtractors.txt), and the
volume of the spheres in a conf-file (conf0 by default).

    python3 solidFraction.py [conf-file]
"""
import math
import sys

WALL_RADIUS = 0.001   # Minkowski radius of the walls (shapes.txt)
XMIN_FACE = 0.0       # inner face of the fixed walls (centre -0.001, radius 0.001)

conf = sys.argv[1] if len(sys.argv) > 1 else "conf0"

# Volume of the spheres: the shape 'Sphere' has a unit radius, so the radius of
# each sphere is its homothety (4th column of its line)
lines = open(conf).read().split("\n")
k = [i for i, l in enumerate(lines) if l.startswith("Particles")][0]
n = int(lines[k].split()[1])
Vs = sum(4.0 / 3.0 * math.pi * float(l.split()[3]) ** 3
         for l in lines[k + 1:k + 2 + n] if l.startswith("Sphere"))


def read(fname):
    return [list(map(float, l.split())) for l in open(fname) if l.strip()]


X, Y, Z = read("xmax.txt"), read("ymax.txt"), read("zmax.txt")
print("# time  Lx  Ly  Lz  solid_fraction")
for x, y, z in zip(X, Y, Z):
    # column 1 is the time, columns 2-4 the position of the wall
    Lx = x[1] - WALL_RADIUS - XMIN_FACE
    Ly = y[2] - WALL_RADIUS - XMIN_FACE
    Lz = z[3] - WALL_RADIUS - XMIN_FACE
    print(f"{x[0]:.4f} {Lx:.5f} {Ly:.5f} {Lz:.5f} {Vs / (Lx * Ly * Lz):.4f}")
