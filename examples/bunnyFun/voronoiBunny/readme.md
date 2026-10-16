## 19. bunnyFun 5: a Voronoi bunny

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

> **This tutorial needs an installation of [Neper](https://neper.info)**, a
> program that is not part of Rockable: see [What is needed](#what-is-needed).

### Goal

The shot of [`bunnyShot`](../bunnyShot) again, on a bunny made of 523 Voronoi
cells instead of cubes. The cells come from **Neper**, a program that
tessellates a volume, and are converted into shapes by **`tess2shape`**. Neper
cannot tessellate the bunny itself: it tessellates a box around it, and a
script keeps the cells whose centre is inside the bunny.

### Notions

- **Neper**: a centroidal Voronoi tessellation of a box (cells of nearly the
  same size).
- The converter **`tess2shape`**: the cells of a tessellation as shapes and
  particles ([converters](https://richefeu.github.io/rockable/converters.html)).
- Keeping the cells inside a shape (`keepInside.sh`, sh + awk), the test of
  `generatePacking:grid_inside` applied to a list of particles.
- A conf-file assembled by a script from several parts (`makeInput.sh`).
- The **time step** of bonded cells, and how to check it.

### What is needed

**Neper** (<https://neper.info>, version 5 here), which is not part of
Rockable. It is built from its sources, with the libraries `gsl` and `nlopt`:

```sh
git clone --depth 1 https://github.com/neperfepx/neper.git
mkdir neper-build && cd neper-build
cmake ../neper/src -DCMAKE_BUILD_TYPE=Release
make -j 8              # the executable neper is in this folder
```

With GCC 16 (`-DCMAKE_C_COMPILER=gcc-16`), add
`-DCMAKE_C_FLAGS="-std=gnu17 -Wno-incompatible-pointer-types"`: GCC 16
compiles C in C23 by default, which a library inside Neper does not accept.
`tess2shape` is built with Rockable.

### Files

| File                 | Content                                                               |
|----------------------|-----------------------------------------------------------------------|
| `../bunny.stl`       | The Stanford bunny, 309 triangles (see [`../readme.md`](../readme.md)) |
| `command.txt`        | The command file of `tess2shape`                                      |
| `keepInside.sh`      | Keeps the cells whose centre is inside the bunny (sh + awk)           |
| `genShapes.txt`      | Commands of `generator` that write `others.shp`: the plate and the ball |
| `params.txt`         | The beginning of the conf-file: the parameters (those of `bunnyShot`) |
| `makeInput.sh`       | Writes `input.txt` and `shapes.txt` from the parts above              |
| `dataExtractors.txt` | `TrackDamage` (`damage.txt`)                                          |
| `fragments.sh`, `energies.sh` | The fragments and the energy balance (as in `bunnyShot`)     |

### Run it

```sh
stl2shape -r 0.001 -c -n Bunny -o bunny.shp ../bunny.stl      # the cutter
neper -T -n 2500 -morpho centroidal -domain "cube(0.21,0.12,0.18)" \
      -transform "translate(-0.114,-0.052,-0.081)" -id 1 -o box   # 70 s, box.tess
tess2shape command.txt                                        # 75 s, cells.txt and cells.shp
sh keepInside.sh bunny.shp Bunny 2 cells.txt cells.shp kept.txt kept.shp
generator genShapes.txt                                       # writes others.shp
sh makeInput.sh                                               # writes input.txt and shapes.txt
rockable input.txt -j 4                                       # 4 minutes (1.5 with -j 8), conf0 to conf40
see
sh fragments.sh conf40
sh energies.sh
```

### What happens

**The tessellation.** `neper -T` tessellates a box of 21 × 12 × 18 cm, the
box that holds the bunny of homothety 2, into 2500 cells. Neper puts its box
at the origin, and `-transform "translate(...)"` moves it around the bunny.
`-morpho centroidal` makes a *centroidal* tessellation: each seed is moved to
the centre of its cell, and the cells end up of nearly the same size (from
1.6 to 2.0 cm³ in the bunny). `-id 1` fixes the random draw: the same
command gives the same cells.

**The conversion.** `tess2shape` reads `box.tess` and the keywords of
`command.txt`, and writes one shape per cell (`cells.shp`) and one particle
per cell (`cells.txt`), in the group 1 and the cluster 0. Each cell is shrunk
by the Minkowski radius (0.2 mm): its faces move inwards by the radius, so
that two neighbouring cells touch without overlapping. The mass properties of
each cell are computed by a Monte-Carlo integration (`MCnstep 5000` points):
most of the 75 seconds.

**The cut.** Neper cannot tessellate the bunny itself: a domain read from a
file (`-domain "file(bunny.obj)"`) makes Neper 5.0.1 crash, even for a cube.
`keepInside.sh` reads the faces of the bunny in `bunny.shp`, scales them by 2,
and keeps the cells whose centre is inside, with the test of
`generatePacking:grid_inside` (rays cast from the centre cross the faces an
odd number of times):

```
523 cells kept out of 2500
```

It writes their lines in `kept.txt` and their shapes in `kept.shp`. The
bunny of `bunny.shp` is not pre-computed: it is only a cutter, in the frame of
the mesh, which is also the frame of the box. `makeInput.sh` then puts
together `params.txt`, the plate, the ball and the cells, lifted by 5.15 cm.

**The glue.** `stickBCM 1e-4` finds the faces of the neighbouring cells, and
glues them: 2781 interfaces, with 5.2 bonds each on average, and 0.152 m² in
all. All the cells are glued in
one piece: the centroidal cells are large enough to keep the ears.

**The shot.** The same ball, at the same place and the same speed as in
`bunnyShot`:

```
sh fragments.sh conf40
cells 523, fragments 98, single cells 87
the largest fragments (cells): 400 6 6 5 4 3

sh energies.sh
kinetic energy of the ball before the impact   198.486 J
energy to break all the interfaces              15.207 J
energy spent in the cracks (conf40)              6.316 J  (3.2% of the kinetic energy)
kinetic energy left in the ball (conf40)         5.864 J
kinetic energy of the cells (conf40)            74.591 J
```

The damage is 0.42 (0.11 with the cubes), and the body keeps 76% of the
cells (89% with the cubes). The two bunnies are not made of the same
material: the stiffness and the strength of the interfaces are given **per
bond** (`knInnerBond`, `gcInnerBond` per unit area), and the interfaces of the
Voronoi cells have more bonds, of other areas, between cells of another size.
Comparing two tessellations needs the parameters of the bonds to be scaled
with them, as in the DCEM paper.

**The time step.** With bonded cells, the time step must resolve the
vibrations of the lightest cell on its bonds. It is checked by running again
with a smaller one:

```
dt       damage   largest fragment
4e-6     1        1        (the bunny explodes: unstable)
2e-6     0.59     74
1e-6     0.42     400
5e-7     0.42     390
```

`1e-6` gives nearly the same result as `5e-7`: it is the largest time step
that can be trusted here. With a random (not centroidal) tessellation, the
smallest cells are seven times lighter than the median one, and the time step
must be even smaller.

### What to look at

- In `see`, the cells of the head flying away, and the faces of the cracks in
  the body.
- `fragments.sh` on the conf-files during the impact, as in `bunnyShot`.
- The number of cells against the number of seeds: `grep -c '^<' kept.shp`.

### Exercises

1. **A random tessellation.** Remove `-morpho centroidal`: the cells have very
   different sizes (from 0.25 to 5.4 cm³), 555 of them are kept, and a whole
   ear (8 cells) is not glued to the head: the thin parts lose their contact
   with the rest when the cells are cut. Find the time step that gives the same
   result as half of it.
2. **Finer cells.** `-n 5000` in the command of Neper: twice as many cells in
   the bunny, and a better shape. Expect the conversion and the computation
   to take at least twice as long, and check the time step again: the cells
   are lighter.
3. **Shoot the body**, at `0 0.05 0.2` in `makeInput.sh`, and compare the
   pieces with those of `bunnyShot`, exercise 2.
4. **Other shapes.** `keepInside.sh` works with any closed shape: cut the box
   with another mesh converted by `stl2shape` (the rock of
   [`complexShape`](../../complexShape), with a homothety of `1e-3` to turn its
   millimetres into metres), after changing the box of Neper so that it holds
   the new shape.

### Going further

- [`bunnyShot`](../bunnyShot): the same shot on cubic cells.
- [`softBallImpactVoronoi`](../../softParticles/softBallImpactVoronoi): a ball
  of 200 soft Voronoi cells, the case of the DCEM paper.
- The commands of Neper: <https://neper.info/doc/>.
