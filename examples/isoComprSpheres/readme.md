## 3. Preparing a sample: isotropic compression

### Goal

125 spheres, placed on a loose grid in a box of six walls, are compressed
until the pressure on the walls is 1000 Pa. This is the usual first step of a
granular test: building a dense sample at a given confining pressure, which
then serves as the initial state of the real test.

### Notions

- A box of six **walls**, the first six bodies (`nDriven 6`).
- A **servo-controller**, `tritriIsostaticCompression`, that imposes a stress
  on the walls ([servoControllers](https://richefeu.github.io/rockable/servoControllers.html)).
- The application **`generator`**, which writes the shapes and the particles.
- The pre-processing command **`homothetyRange`**, which gives the spheres
  random sizes ([preProcessing](https://richefeu.github.io/rockable/preProcessing.html)).
- **Data extractors** to follow the walls, and the **solid fraction**.
- A conf-file used as the **initial state** of another computation.

### Files

| File                 | Content                                                             |
|----------------------|---------------------------------------------------------------------|
| `input.txt`          | The conf-file: 6 walls and 125 spheres                              |
| `shapes.txt`         | The shapes `Sphere` (unit radius) and `x-wall`, `y-wall`, `z-wall`   |
| `drivingSystem.txt`  | The servo that compresses the box                                   |
| `dataExtractors.txt` | Records the positions of the three moving walls                     |
| `genShapes.txt`      | Commands of `generator` that write `shapes.txt`                      |
| `genSample.txt`      | Commands of `generator` that write the walls and the spheres         |
| `solidFraction.py`   | Computes the solid fraction from the recorded wall positions         |

### Run it

```sh
rockable input.txt           # about 20 seconds
see conf100
python3 solidFraction.py     # the solid fraction during the compression
```

### What happens

**Building the system.** The shapes and the particles were written by
`generator`:

```sh
generator genShapes.txt > shapes.txt
generator genSample.txt > particles.txt
```

`genShapes.txt` makes a sphere of unit radius and the three wall shapes of a
box of 5 cm. `genSample.txt` places the six walls (`generatePacking:wallBox`),
then 5 × 5 × 5 spheres on a grid (`generatePacking:grid`). The lines of
`particles.txt` were pasted after the keyword `Particles` of `input.txt`.

**Random sizes.** The spheres are all written with the same homothety. The
last line of `input.txt`,

```
homothetyRange 6 130   2.0e-3 4.0e-3  0
```

gives the particles 6 to 130 a homothety drawn between 2 and 4 mm; since the
shape `Sphere` has a unit radius, this is their radius. The last `0` makes the
draw reproducible (`1` would seed it on the clock).

**The parameters.** There is no gravity (`gravity 0 0 0`): the sample is only
loaded by the walls. Group 0 is the walls, group 1 the spheres, with a
friction coefficient of 0.2.

**The servo.** `drivingSystem.txt` holds

```
Servo tritriIsostaticCompression 0 1 2 3 4 5
1000
```

The six numbers are the walls `Xmin Xmax Ymin Ymax Zmin Zmax`. The three `min`
walls stay fixed, and a force is applied on the three `max` walls, recomputed
at every step as the pressure (1000 Pa) times the current area of each face.
The walls move inwards, the spheres get in contact, and the box stops
shrinking when the contacts carry the pressure.

The box goes from 50 mm to about 31.7 × 32.5 × 30.0 mm in about 0.1 s, then
hardly moves. The solid fraction (the volume of the spheres over the volume of
the box) goes from 0.13 to **0.52**: a rather loose packing, because of the
friction and of the fast compression.

### What to look at

- In `see`, compare `conf0` and `conf100`. The key `d` hides the walls, and
  the key `2` colours the spheres by their velocity.
- `python3 solidFraction.py` prints the size of the box and the solid fraction
  at each record: it reaches its final value around `t = 0.1 s`.
- `gnuplot checkplots.txt`: the kinetic energy drops by several orders of
  magnitude once the sample is compressed.

### Exercises

1. Set the friction between spheres to zero (`muContact 1 1 0.0`): the sample
   is denser (solid fraction about 0.54).
2. Compress at 10 000 Pa instead of 1000 Pa. How much does the solid fraction
   change?
3. Make the sizes more dispersed: `homothetyRange 6 130 1.0e-3 4.0e-3 0`.
4. Use the compressed sample for a biaxial test. Copy `conf100` to a new
   folder, with `shapes.txt` and `dataExtractors.txt`, set `tmax 1.5` in it,
   and write this `drivingSystem.txt`:

   ```
   Servo tritriBiaxialCompression 0 1 2 3 4 5
   1000 0.01
   ```

   then run `rockable conf100`. The top wall now comes down at 0.01 m/s, while
   the lateral pressure stays at 1000 Pa: the sample shortens along `y`,
   widens along `x` and `z`, and its solid fraction decreases (it dilates).

### Going further

- [`prepro_cpptools`](../prepro_cpptools): more about `generator`.
- [`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic): the same compression
  in a periodic cell, without walls.
- The other servos:
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html).
