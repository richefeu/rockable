## 5. A complete system built with `generator`: dominoes and a glued wall

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

Six dominoes, each 1.5 times larger than the previous one, topple one after the
other, and the last one, of nearly 5 kg, falls on a wall of 125 small bricks
glued together. The whole system, shapes and input file, is written by a single
command file of `generator`.

### Notions

- **`generator`** to write a whole system: `open`, `print`, `close`,
  `generateShape:cuboid`, `generateShape:rectangle_xz`,
  `generatePacking:grid_clust`
  ([generator](https://richefeu.github.io/rockable/generator.html)).
- **Clusters** of glued particles (`stickVerticesInClusters`) and the law
  `StickedLinks`, as in [`stickVerticesInClusters`](../stickVerticesInClusters).
- The **time step**, and the numerical stability of the damping of the bonds.

### Files

| File           | Content                                                               |
|----------------|-----------------------------------------------------------------------|
| `generate.txt` | Command file of `generator`: writes `shapes.txt` and `input.txt`       |

### Run it

```sh
generator generate.txt       # writes shapes.txt and input.txt
rockable input.txt -j 4      # about 5 minutes
see
```

### What happens

**The shapes.** The first block of `generate.txt`, between `open shapes.txt`
and `close`, writes the shapes:

```
generateShape:rectangle_xz SOL0 1e-2 100e-2 30e-2
generateShape:cuboid Domino_1   0.00125  0.005  0.04  0.02
generateShape:cuboid Domino_2   0.001875 0.0075 0.06  0.03
...
generateShape:cuboid Domino_6   0.009492188 0.03796875 0.30375 0.151875
```

a floor of 1 m × 30 cm (radius 1 cm), and six cuboids (radius, then the three
sizes). Each domino is 1.5 times larger than the previous one: the smallest is
4 cm high and weighs 11 g, the sixth is 30 cm high and weighs 4.7 kg. A domino
can topple another one about 1.5 times larger, so that a small push ends up
moving a large mass: the lines in comment go up to `Domino_19`, 59 m high.

**The input file.** The second block, `open input.txt`, writes the input file
with `print` lines: the parameters, then the particles. The floor is the only
driven body. Each domino has its own cluster number (`1` to `6`), so that they
are not glued together. The first one is given an angular velocity of
-10 rad/s about `z`: it is the push.

**The wall.** The command

```
generatePacking:grid_clust Domino_1  0.55 0.0 -0.05   0.025 0.2 0.1   5 5 5    1  8  1.0  0
```

places 5 × 5 × 5 bricks of the shape `Domino_1` on a grid starting at
`(0.55, 0, -0.05)` and of size `0.025 × 0.2 × 0.1`: the bricks touch each
other. They are in group 1, all in cluster 8, with a homothety of 1 and no
random orientation (last `0`). The line `stickVerticesInClusters 0.001`, at
the end of the input file, glues the vertices of the bricks of a same cluster
closer than 1 mm: 300 interfaces of 4 bonds. The bonds break at 5 N
(`fn0InnerBond`, `ft0InnerBond`).

**The dynamics.** The dominoes fall one after the other, the larger ones
taking longer: the sixth one starts to tip at about `t = 1.8 s`, and hits the
wall at `t = 1.97 s`. The wall is shattered: all its interfaces break
within 0.02 s, and the bricks are thrown along `x` at about 2.4 m/s.

**The time step and the damping of the bonds.** The bricks are thin and light,
and each one is held by up to 24 bonds of stiffness `1e7`. With
`en2InnerBond 1 1 0.02` (a strong damping), the explicit integration of the
viscous force of the bonds is unstable: the rotation of the bricks grows from
round-off, and the wall falls apart by itself within a millisecond, before any
impact. With `en2InnerBond 1 1 0.5`, the wall stays intact, even with a time
step of `2e-6`. See the warning in
[forceLaws](https://richefeu.github.io/rockable/forceLaws.html) (`StickedLinks`).

### What to look at

- In `see`, the chain of dominoes, then the impact on the wall.
- The number of interfaces still glued, in each conf-file:

  ```sh
  for i in $(seq 90 105); do echo "conf$i $(grep '^Interfaces' conf$i)"; done
  ```
- The kinetic energy, which grows at each domino: `kineticEnergy.txt`.

### Exercises

1. Make the glue ten times stronger (`fn0InnerBond 1 1 50.` and
   `ft0InnerBond 1 1 50.0`, in `generate.txt`), run `generator` again, and
   compare what remains of the wall.
2. Set `print en2InnerBond 1 1 0.02` and `print tmax 0.01`, and count the
   interfaces in the conf-files: the wall falls apart without being touched.
   Then set `print dt 5e-7` as well: it holds again.
3. Uncomment `Domino_7` (the shape and the particle), and change
   `print Particles 132` into `print Particles 133`: a seventh domino, 46 cm
   high, stands between the sixth one and the wall.
4. Give random orientations to the bricks (last `0` of `grid_clust` set to
   `1`): they now overlap each other, only 5 interfaces are made, and the
   overlaps throw the bricks apart at about 10 m/s within a millisecond. A
   grid is only a valid packing when its cells are larger than the particles.

### Going further

- [`prepro_cpptools`](../prepro_cpptools): `generator` step by step.
- [`stickVerticesInClusters`](../stickVerticesInClusters): the bonds and
  their failure criterion.
