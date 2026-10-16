## 20. bunnyFun 6: a bunny in blocks

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

The bunny is cut into 85 blocks by **`stl2blocks`**: vertical slices every
3 cm, each slice cut into bars every 3 cm, and each bar cut in three by a
tilted plane. The blocks are glued, and the bunny stands on its feet on a
plate. A steel cannonball, 12 cm wide, rolls towards it at 3 m/s, and breaks
it into 25 fragments.

### Notions

- The converter **`stl2blocks`**: a mesh cut into blocks along a **tree of
  plane cuts**, each block with its own direction and its own spacing of the
  planes ([converters](https://richefeu.github.io/rockable/converters.html)).
- The Minkowski radius in a cut: the faces of the mesh moved inwards by `R`,
  and a slice of `2R` removed by each plane, so that the rounded blocks fill
  the mesh and touch each other.
- The **vertices for the glue**, added by `stl2blocks` on both sides of each
  plane, so that `stickBCM` finds vertices that face each other.
- A **rolling** ball: a velocity and an angular velocity such that `v = ω R`.
- The fragments measured by their **volume** (`fragments.sh`, sh + awk): the
  blocks do not have the same size.

### Files

| File                 | Content                                                               |
|----------------------|-----------------------------------------------------------------------|
| `../bunny.stl`       | The Stanford bunny, 309 triangles (see [`../readme.md`](../readme.md)) |
| `cut.txt`            | The command file of `stl2blocks`: the tree of cuts                    |
| `genShapes.txt`      | Commands of `generator` that write `others.shp`: the plate and the cannonball |
| `params.txt`         | The beginning of the conf-file: the parameters (those of `bunnyShot`, with a weaker glue) |
| `makeInput.sh`       | Writes `input.txt` and `shapes.txt`; the speed of the cannonball is its argument |
| `dataExtractors.txt` | `TrackDamage` (`damage.txt`)                                          |
| `fragments.sh`       | The fragments of a conf-file, by volume                               |
| `energies.sh`        | The energy balance of the impact                                      |

### Run it

```sh
stl2blocks cut.txt          # 2 seconds: blocks.shp and blocks.txt
generator genShapes.txt     # writes others.shp
sh makeInput.sh             # writes input.txt and shapes.txt (cannonball at 3 m/s)
rockable input.txt -j 4     # 2 minutes, conf0 to conf60
see
sh fragments.sh
sh energies.sh
```

### What happens

**The cuts.** In `cut.txt`, the bunny is read with a homothety of 2 (21 cm
long), lifted so that it stands on the plate, and given a Minkowski radius of
1 mm. Three lines of cuts make the tree:

```
cut root      1 0 0      spacing 0.03 0
cut root.*    0 0 1      spacing 0.03 0.015
cut root.*.*  0 1 0.3    number 3
```

The whole bunny (`root`) is cut by the planes `x = 0.03 k`, into slices
numbered `0`, `1`, ... along `x`. Each slice (`root.*`: `*` is any slice) is
cut along `z`, every 3 cm; and each bar (`root.*.*`) is cut into three blocks
of the same thickness, along a direction tilted by 17° from the vertical.
`stl2blocks` prints the tree:

```
root: direction 1 0 0, 7 planes, slabs 0 1 2 3 4 5 6 7
  root.0: direction 0 0 1, 2 planes, slabs 0 1 2
    root.0.0: direction 0 0.957826 0.287348, 2 planes, slabs 0 1 2
    ...
```

then one line per block, with its path (`root.0.0.1` is the block 1 of the
bar 0 of the slice 0) and its volume. The faces of the mesh are first moved
inwards by 1 mm, and each plane removes a slice of 2 mm: rounded by 1 mm, the
blocks fill the bunny, and two neighbours touch. At the surface, a plane can
cut a sliver: the parts smaller than `minVolume` (0.3 cm³) are removed.

```
  20 tiny parts removed (core volume <= minVolume)
> 1717 vertices added in the faces of the cuts, for the glue (stickBCM)
  area of the faces of the cuts: 0.0776607
  smallest block: block_33, volume 7.41163e-07 (it sets the time step)
> 85 blocks, total volume 0.000941807 (98.6888% of the mesh)
```

The mesh is repaired first: the bunny has one triangle stuck to the others by
a single edge (a fin), which is removed.

**The glue.** `stickBCM` bonds the vertices of two blocks that face each
other. Two neighbours are cut by different planes after the plane between
them, and their corners do not face each other; `stl2blocks` adds to both the
corners of the part of the plane they share (1717 vertices here, in the faces,
which do not change the shapes):

```
stickBCM: 233 interfaces created, with 10.5966 bonds per interface on average, total area 0.0719533
```

**Standing.** The lowest block starts 0.1 mm above the plate: the bunny lands
and stands still (its blocks move at less than 1 mm/s) until the cannonball
arrives.

**The cannonball.** `makeInput.sh` writes its line by hand: a sphere of 6 cm
(7.1 kg) at 4 cm from the bunny, with the velocity `0 0 -3` and the angular
velocity `-50 0 0`. Its contact point with the plate does not slide
(`v = ω R`): it rolls, and keeps rolling at 3 m/s. It touches the bunny at
`t = 15 ms`.

**The impact.** Within 25 ms, the bunny breaks; the cannonball goes through,
and rolls over the debris:

```
sh fragments.sh
blocks 85, fragments 25, single blocks 16
the largest fragments (% of the volume, blocks):  49.4 (32)  13.3 (11)  7.5 (9)  6.3 (3)  3.9 (2)

sh energies.sh
kinetic energy of the cannonball before the impact   44.461 J
energy to break all the interfaces                    2.159 J
energy spent in the cracks (conf60)                   0.970 J  (2.2% of the kinetic energy)
kinetic energy left in the cannonball (conf60)       22.141 J
kinetic energy of the blocks (conf60)                 4.745 J
```

The damage (the broken fraction of the area of the interfaces, `damage.txt`)
is 0.45. The kinetic energy of the cannonball counts its rotation (12.7 J of
the 44.5 J). The glue is weaker than in `bunnyShot` (`gcInnerBond 1 1 15`
instead of 50): breaking all the interfaces would take 2.2 J, and the cracks
take 1 J.

**The time step.** The smallest block (0.74 cm³, 2 g) sets it, on its bonds:
`stl2blocks` names it. `dt 1e-6` gives nearly the same result as `5e-7` (a
damage of 0.45 against 0.42, a largest fragment of 49% against 48% of the
volume).

### What to look at

- In `see`, the key `3` gives each shape its own colour: the blocks, and the
  slices, become visible. Look at the bars of the slices, tilted, in the
  fragments.
- `fragments.sh` on the conf-files during the impact (`conf8`, `conf10`,
  `conf15`, `conf20`): the bunny comes apart between 15 and 40 ms.
- The damage against time:

  ```sh
  gnuplot -p -e "plot 'damage.txt' using 1:2 with lines title 'damage'"
  ```

### Exercises

1. **The speed of the cannonball.** `sh makeInput.sh 2`, then `rockable`: at
   2 m/s, the cannonball has 19.8 J instead of 44.5. The damage is 0.18, and
   the bunny keeps 81% of its volume in one piece; the cannonball climbs onto
   it (its centre ends 3.5 cm higher).
2. **A tougher glue.** `gcInnerBond 1 1 50` in `params.txt` (that of
   `bunnyShot`), at 3 m/s: the damage is only 0.09, and the bunny keeps 95%
   of its volume in one piece.
3. **Without the vertices for the glue.** Add `gluePoints 0` to `cut.txt`:
   `stickBCM` finds 124 interfaces, with 2.7 bonds each. An interface with
   fewer than three bonds has no area, and breaks at the first step; two
   blocks without facing vertices have no interface at all. Before the
   cannonball arrives, the bunny is already in 28 pieces, held by friction
   (yet the damage is zero: these interfaces had no area); the cannonball
   scatters it into 54 fragments.
4. **Smaller blocks.** Cut every 2 cm (`spacing 0.02 0` and
   `spacing 0.02 0.01`): 142 blocks, 97% of the volume (more slivers are
   removed). The damage is 0.24, and the bunny breaks into two large pieces
   (57% and 25% of the volume).
5. **A cannonball that slides.** In `makeInput.sh`, write `0` instead of
   `-V / 0.06` (its angular velocity): it slides on the plate instead of
   rolling, with 31.8 J instead of 44.5. The damage is 0.39 (0.45 when it
   rolls), and the largest fragment keeps 57% of the volume: the energy of
   rotation does little to break the bunny.
6. **Another tree.** Add `cut root.3.* 1 1 1 random 2 7` to `cut.txt`: the
   bars of the slice 3 are cut by two planes at random positions, normal to
   the diagonal `1 1 1`, instead of the tilted planes (when several lines
   match a block, the one with the fewest `*` is used). The bunny has 84
   blocks; the damage is 0.42, and the largest fragment keeps 48% of the
   volume.

### Going further

- [`bunnyShot`](../bunnyShot) and [`voronoiBunny`](../voronoiBunny): the
  bunny cut into cubes and into Voronoi cells.
- [`fragmentation`](../../fragmentation): the law `BCM` and its parameters.
- `stl2blocks` and its keywords:
  [converters](https://richefeu.github.io/rockable/converters.html).
