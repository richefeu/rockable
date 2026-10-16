`stl2blocks` cuts a closed STL mesh (binary or ascii) into blocks, along families of parallel planes, and writes the blocks as Rockable shapes and particles, ready to be glued by `stickBCM`.

```sh
stl2blocks command.txt
```

The cuts form a tree: the whole mesh (`root`) is cut along a direction into slabs (`root.0`, `root.1`...), and each slab can be cut again along its own direction, with its own spacing of the planes. The Minkowski radius R is given first: the faces of the mesh are moved inwards by R, and each cutting plane removes a slice of thickness 2R, so that the blocks, rounded by R, fill the mesh and touch each other.

```
stlFile          bunny.stl
homothety        2
MinkowskiRadius  0.001
cut root      1 0 0      spacing 0.02 0      # planes x = 0 + k 0.02
cut root.*    0 1 0.3    spacing 0.02 0.005  # every slab of root
cut root.*.*  0.2 -0.3 1 number 5            # 5 slabs in each block
cut root.3.*  1 1 1      random 2 3          # but these ones: 2 random planes (seed 3)
shapesFile       blocks.shp
inputFile        blocks.txt
ParticlesCluster 1
```

The planes are given by `number <n>`, `spacing <s> <offset>`, `positions <k> <p1>...<pk>`, `fractions <k> <f1>...<fk>` or `random <k> <seed>`. A `*` in a path matches any slab; the rule with the fewest `*` wins.

The other keywords are `translation`, `shapeName`, `ParticlesGroup`, `minVolume`, `MCnstep` and `gluePoints`. The documentation (page *Converters*) gives the details: how the vertices for the glue are added on both sides of each plane, how the mesh is repaired, and the limits (parts thinner than 2R).

`polyhTool.hpp` is a copy of the one of toofus, in the namespace `stl2blocks`, extended to polyhedra that are not convex (cut by a plane with faces of several loops, with holes, concave; offset of the faces; separation of the parts).
