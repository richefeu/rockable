## bunnyFun: tutorials with the Stanford bunny

<p align="center"><img src="bunnyPile/animation.gif" alt="a pile of bunnies"/></p>

The Stanford bunny is the most famous mesh of computer graphics: a ceramic
rabbit scanned in 1994 by Greg Turk and Marc Levoy at Stanford. These
tutorials use it as a grain, as a container and as a block to break, in order
to put the pre-processing tools of Rockable to work on a shape that is neither
convex nor simple: `stl2shape`, `shapeSurvey`, `shape2mesh`, and the packing
commands of `generator`.

| Folder                   | What it shows                                                                   | Time            |
|--------------------------|---------------------------------------------------------------------------------|-----------------|
| [`bunnyPile`](bunnyPile) | A pile of bunnies: from the mesh to the shape, packing any shape with `generatePacking:OBBs`, a tube lifted | 50 s (4 threads) |
| [`bunnyInBunny`](bunnyInBunny) | Bunnies in a bunny: a shape edited into a shell with a hole (`isSurface`), a driven rotation, counting with a ray | 1.5 min (4 threads) |
| [`bunnyShot`](bunnyShot) | A shot in the head: a bunny cut into glued cubic cells (`generatePacking:grid_clust_inside`, `BCM`), fragments and energy | 45 s (4 threads) |

### The mesh

`bunny.stl` is a binary STL of the bunny, decimated to 309 triangles (157
vertices) with MeshLab, from the 69 451 triangles of the original scan. It is
closed (each edge belongs to two triangles), and it is about 10 cm long, in
metres.

The original data come from the
[Stanford 3D Scanning Repository](https://graphics.stanford.edu/data/3Dscanrep/)
of the Stanford Computer Graphics Laboratory, which allows any use of the
bunny, as long as the source is acknowledged.
