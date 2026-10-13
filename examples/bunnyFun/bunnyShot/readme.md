## 17. bunnyFun 3: a shot in the head

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

The Stanford bunny is cut into 956 cubic cells of 1 cm, glued face to face,
and sits on a plate. A steel ball of 3 cm hits its head at 60 m/s. The head
bursts into small fragments, and the body stays whole. The cells are made by a
**cookie cutter**: a grid of cubes from which `generator` keeps only the cubes
inside the bunny.

### Notions

- The packing commands **`generatePacking:grid_inside`** and
  **`generatePacking:grid_clust_inside`**: a grid cut by any closed shape
  ([generator](https://richefeu.github.io/rockable/generator.html)).
- A particle written by hand in a command file of `generator` (`print`), with
  its velocity, and counted with **`incrementNoParticles`**.
- The bonded cells (`BCM`, `stickBCM`) of [`fragmentation`](../../fragmentation),
  on a shape that is not a box.
- The **fragments**, found by following the interfaces left in a conf-file
  (`fragments.sh`, sh + awk), and the energy balance of the impact.

### Files

| File                 | Content                                                               |
|----------------------|-----------------------------------------------------------------------|
| `../bunny.stl`       | The Stanford bunny, 309 triangles (see [`../readme.md`](../readme.md)) |
| `genShapes.txt`      | Commands of `generator` that write `shapes.txt`: the cell, the plate, the ball |
| `genInput.txt`       | Commands of `generator` that write `input.txt`: the plate, the ball, the cells |
| `dataExtractors.txt` | `TrackDamage`: the broken fraction of the interfaces (`damage.txt`)   |
| `fragments.sh`       | The number and the sizes of the fragments in a conf-file             |
| `energies.sh`        | The energy balance of the shot                                       |

### Run it

```sh
stl2shape -r 0.001 -c -n Bunny -o bunny.shp ../bunny.stl   # the cutter
shapeSurvey bunny.shp -P bunny.shp
generator genShapes.txt                                    # writes shapes.txt
generator genInput.txt                                     # writes input.txt
rockable input.txt -j 4                                    # 45 seconds, conf0 to conf40
shape2mesh shapes.txt --rmsh -s
see
sh fragments.sh conf0
sh fragments.sh conf40
sh energies.sh
```

### What happens

**The cutter.** In `genInput.txt`, the line

```
generatePacking:grid_clust_inside bunny.shp Bunny 2  Cell  -0.12 -0.06 -0.09  0.22 0.14 0.19  22 14 19  1 2 1 0
```

reads the shape `Bunny` in `bunny.shp`, scales it by 2 (21 cm long), and puts
it at the origin, in the frame of the shape file. The rest of the line is that
of `generatePacking:grid_clust`: a grid of 22 × 14 × 19 positions, 1 cm apart,
in the box that holds the scaled bunny. Only the positions inside the bunny
are kept: from each one, rays are cast, and a position is inside when the
rays cross the faces of the bunny an odd number of times.

```
@generatePacking:grid_clust_inside, Number of added bodies: 956 (out of 5852 in the grid)
```

The cells fill 956 cm³, 88% of the volume of the bunny (1090 cm³): the
surface is a staircase, which cuts the cubes on the edge at their centre. All
the cells are in the cluster 2: `stickBCM 1e-4`, at the end of `input.txt`,
glues them face to face, with 2392 interfaces. The thin parts of the bunny do
not all survive: three cells at the tip of an ear touch the rest only by an
edge, and have no interface with it (`fragments.sh conf0` finds two
fragments, of 953 and 3 cells).

The transformation `add_translation 0 0.0505 0` is applied to the positions
kept, after the cut: it lifts the bunny so that its lowest cells lie on the
plate.

**The ball.** `generator` has no command for a single particle with a
velocity: the line of the ball is written by hand, with `print`, and
`incrementNoParticles 1` adds it to the number of particles written by
`Particles:placeHolder`. The ball (group 2, steel, 110 g) starts 20 cm in
front of the head, at 60 m/s: 198 J.

**The shot.** The ball reaches the head at `t = 1.8 ms`. `damage.txt` gives the
broken fraction of the interfaces, 0.113 at the end. `fragments.sh` follows,
in a conf-file, the interfaces that are left: two cells linked by an
interface are in the same fragment.

```
sh fragments.sh conf40
cells 956, fragments 84, single cells 72
the largest fragments (cells): 850 12 3 3 2 2
```

The head is gone: about a hundred cells, most of them single, fly away, and
the body keeps 850 cells. The energy balance:

```
kinetic energy of the ball before the impact   198.486 J
energy to break all the interfaces              15.309 J
energy spent in the cracks (conf40)              1.728 J  (0.9% of the kinetic energy)
kinetic energy left in the ball (conf40)        10.990 J
kinetic energy of the cells (conf40)           168.381 J
```

Less than 1% of the energy of the ball opens cracks: it is mostly given to the
106 cells of the head, which fly away at 16 m/s on average (much faster for a
few of them), while the body slides back at less than 1 m/s.

### What to look at

- In `see`, the head bursting, and the body that hardly moves (the keys `0` to
  `5` change the colours; with `1`, random colours, the cells are easier to
  tell apart). The skins
  (`shapes.rmsh`) give the cells their rounded edges.
- The damage against time:

  ```sh
  gnuplot -p -e "plot 'damage.txt' using 1:2 with lines title 'damage'"
  ```
- `fragments.sh` on the conf-files during the impact (`conf5`, `conf10`...):
  the head comes apart within a few milliseconds.

### Exercises

1. **The size of the bunny.** Only `generator` is needed: change the homothety
   of the cutter (the `2` after `Bunny`), and widen the grid so that it holds
   the bunny (for instance `-0.18 -0.09 -0.14  0.33 0.21 0.28  33 21 28`). The
   cells fill 118, 403, 956 and 3206 cm³ for homotheties of 1, 1.5, 2 and 3:
   always 87 to 88% of the volume of the bunny. With 1 cm cells, a bunny of
   homothety 1 (10 cm) keeps only one cell of its ears.
2. **Shoot the body.** Put the ball at `0 0.05 0.2` in `genInput.txt` (it now
   hits the flank, in the middle of the body). At 30 m/s, it digs a crater and
   knocks a few cells off the far side (damage 0.07, the bunny keeps 911
   cells). At 60 m/s, the bunny breaks into three large pieces, of 272, 231
   and 153 cells (damage 0.40); at 100 m/s, into four, of 146, 143, 115 and 87
   cells (damage 0.49).
3. **A tougher material.** `print gcInnerBond 1 1 200` in `genInput.txt`. The
   head still bursts, into a few larger pieces (21 and 17 cells), and the body
   keeps 828 cells (damage 0.10). The cracks now take 6.2 J, 3.1% of the
   energy of the ball.
4. **A scattered toughness.** After `print stickBCM 1e-4`, add
   `print setVariableStickParams Gc inner 50 4 0`: each interface draws its
   own `Gc` from a Weibull law (see [`fragmentation`](../../fragmentation),
   exercise 7). The weak interfaces break first: 99 fragments, a body of 807
   cells (damage 0.13).
5. **The ball.** Make it heavier (a radius of 2 cm in `genShapes.txt`), or
   aim at the ears.

### Going further

- [`fragmentation`](../../fragmentation): the law `BCM` in detail, and the
  post-processor `ClusterGranulo`.
- [`bunnyPile`](../bunnyPile) and [`bunnyInBunny`](../bunnyInBunny): the
  bunny as a grain and as a container.
- All the packing commands of `generator`:
  [generator](https://richefeu.github.io/rockable/generator.html).
