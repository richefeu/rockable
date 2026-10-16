## 12. A rockfall: three blocks down a slope

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

Three blocks of the same volume (about 1 m³) but of different shapes, a
sphere, a cube and a slab, are released side by side on a slope of 40°, 15 m
from its top, and stop on the horizontal run-out zone below. How far each one
goes depends on its shape: this is the kind of question a rockfall study
asks, and what Rockable was first written for.

### Notions

- A **slope**: a driven plate turned by a quaternion.
- The force law **`Avalanche`**: hysteretic normal force, friction, and
  rolling resistance ([forceLaws](https://richefeu.github.io/rockable/forceLaws.html)).
- The **run-out** of a block, measured with `TrackBody`, and the data
  extractor **`TrackRockfall`**, which stops the computation when the block
  stops ([dataExtractors](https://richefeu.github.io/rockable/dataExtractors.html)).

### Files

| File                 | Content                                                        |
|----------------------|----------------------------------------------------------------|
| `input.txt`          | The conf-file: the slope, the run-out zone and the three blocks |
| `genShapes.txt`      | Commands of `generator` that write `shapes.txt`                |
| `shapes.txt`         | The shapes `Slope`, `Floor`, `Sphere`, `Cube` and `Slab`        |
| `dataExtractors.txt` | Follows the three blocks (`sphere.txt`, `cube.txt`, `slab.txt`) |
| `runout.sh`          | Maximum speed, stop time and run-out of each block             |

### Run it

```sh
rockable input.txt     # less than a second
see
sh runout.sh
```

### What happens

**The geometry.** `genShapes.txt` makes two plates (`generateShape:rectangle_xz`),
40 × 12 m and 60 × 12 m with a radius of 0.1 m, a sphere of radius 0.62 m, a
rounded cube of 1 m and a slab of 2 × 0.5 × 1 m: the three blocks have nearly
the same volume. The plates are the two driven bodies. The slope is centred
on the origin and turned by -40° about `z`, which is the quaternion

```
Q = (cos(-20°), 0, 0, sin(-20°)) = (0.939693, 0, 0, -0.342020)
```

(a rotation of angle θ about the unit axis **u** is the quaternion
`(cos θ/2, sin θ/2 u)`). Its foot is at `(20 cos 40°, -20 sin 40°)`, that is
`(15.32, -12.86)`, where the run-out zone starts. The blocks are placed 15 m
up the slope, at a small distance above it, with the same rotation as the
slope for the cube and the slab, so that they lie flat on it.

**The law.** `Avalanche` has no viscous damping: the normal force is
hysteretic, stiff when the bodies get closer and softer when they separate,
so that a fraction `en2Contact` of the energy is given back (here 0.1). The
tangential force is capped by the friction (`muContact 0.5`), and the moment
by `murContact` (here 0.3, dimensionless in this law) times the lever arm and
the normal force: it is a rolling resistance, without which a sphere would
roll for ever.

**What the blocks do.**

```
# block  max speed (m/s)  stop time (s)  run-out beyond the foot (m)
sphere      13.99            9.6              24.03
cube        13.06            7.0               6.92
slab        13.02            7.1               0.19
```

The cube and the slab **slide** down the slope, without rotating, since
`tan 40° = 0.84` is larger than the friction 0.5: their acceleration is
`g (sin 40° - 0.5 cos 40°) = 2.55 m/s²`, and after 35 m they reach
`13.4 m/s` (13.1 in the computation). On the flat zone, the friction stops
them: the slab at once, the cube after tipping over once. The sphere **rolls**,
and the rolling resistance takes it 24 m beyond the foot.

### What to look at

- In `see`, the three blocks, from the side (the key `d` hides the driven
  plates).
- The trajectories, from the files of `TrackBody` (time, position, velocity,
  orientation, angular velocity...: see `extractedDataDoc.txt`):

  ```sh
  gnuplot -p -e "plot 'sphere.txt' using 2:3 with lines title 'sphere', 'cube.txt' using 2:3 with lines title 'cube', 'slab.txt' using 2:3 with lines title 'slab'"
  ```
- The speed against time (`using 1:(sqrt($5**2+$6**2+$7**2))` in gnuplot).

### Exercises

1. Remove the rolling resistance (`murContact 0 1 0.0`): the sphere rolls off
   the end of the run-out zone and falls. The cube and the slab go further
   too (13 m and 1 m): the resistant moment also acts on them.
2. Lower the friction (`muContact 0 1 0.3`): the slab now slides 8 m beyond
   the foot, and the cube 15 m.
3. Follow only the sphere with `TrackRockfall` (`dataExtractors.txt`):

   ```
   TrackRockfall 2 0.05 0.1 sphere.txt 1000
   ```

   The computation stops by itself when the speed of the sphere is below
   0.05 m/s and its angular velocity below 0.1 rad/s (at `t = 9.6 s`), or when
   it leaves the box of the driven bodies: with exercise 1, at `t = 8 s`.
4. Make the slope steeper or gentler: the angle appears in the quaternion of
   the slope and of the cube and the slab, in the position of the run-out zone
   and in the positions of the blocks. Below 26.6° (`tan α = 0.5`), the cube
   and the slab do not slide at all.
5. Release the blocks 3 m higher, with `en2Contact 0 1 0.5`: they bounce when
   they land, but end nearly at the same place.

### Going further

- A real slope is rough: a few driven obstacles, or a mesh of the ground
  converted with `stl2shape` (see [`complexShape`](../complexShape)), make the
  blocks bounce and turn.
- The other force laws, and the dissipation:
  [forceLaws](https://richefeu.github.io/rockable/forceLaws.html),
  [dissipation](https://richefeu.github.io/rockable/dissipation.html).
