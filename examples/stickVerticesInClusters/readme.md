## How to build a breakable cluster?

This example was set as a benchmark for PhD Thesis of Marta Stasiak.


First, in the input file, you need to use the force law named `StickedLinks`

```
forceLaw StickedLinks
```

The required parameters can be of 3 types : simple contacts, `Inner` glue or `Outer` glue. They are set by specifying a pair of group-id.

```
knContact 0 1 0.8e7
ktContact 0 1 0.8e7
en2Contact 0 1 0.1
muContact 0 1 0.5
krContact 0 1 0.0
murContact 0 1 0.0

knContact 0 0 0.8e7
ktContact 0 0 0.8e7
en2Contact 0 0 0.1
muContact 0 0 0.5
krContact 0 0 0.0
murContact 0 0 0.0

knInnerBond 0 0 0.8e7
en2InnerBond 0 0 0.1
ktInnerBond 0 0 0.8e7
fn0InnerBond 0 0 85.
ft0InnerBond 0 0 170.0
powInnerBond 0 0 2.0
```

When defining the particles, the group-id but also the cluster-id need to be set. for example, in the following line, the group-id is `0` (value just after the shape name) and the cluster-id is `0` (second value after the shape name). 

```
ShellSector8_N 0 0 1 ...
```

Finally, at the end of the input file, the preprocessing command `stickVerticesInClusters` is employed to glue (or to stick) the close vertices that belong to the same cluster (i.e., they share the same cluster-id). The value given just after the keyword is the distance below which the glue is activated.

```
stickVerticesInClusters 0.0004
```

It is important to notice that, according to the parameters used for gluing the vertices, the solid link created may be directly broken. But actually the relative displacements are set to zero so that all forces are also null.

One other thing to notice, is that the glue is activated only in-between the vertices, and the normal direction is defined from the center-to-center branch (the surface normal is _**not**_ explicitly used).