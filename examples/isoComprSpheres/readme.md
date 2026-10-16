## 3. Preparing a sample: isotropic compression

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

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
| `solidFraction.sh`   | Computes the solid fraction from the recorded wall positions         |
| `state.txt`          | Command file of `postpro`: the same, and more, from the conf-files   |

### Run it

```sh
rockable input.txt           # about 20 seconds
see conf100
sh solidFraction.sh          # the solid fraction during the compression
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
- `sh solidFraction.sh` prints the size of the box and the solid fraction at
  each record: it reaches its final value around `t = 0.1 s`. The script only
  uses `awk`; read it to see how the files of the data extractors are combined.
  To plot the solid fraction:

  ```sh
  sh solidFraction.sh > phi.txt
  gnuplot -p -e "plot 'phi.txt' using 1:5 with lines title 'solid fraction'"
  ```
- The same with `postpro`, from the conf-files: `postpro state.txt` runs the
  post-processor `SampleState` (see [`postProcessing`](../postProcessing)) on
  `conf0` to `conf100`, and writes `SampleState.txt`: the sizes of the box
  (columns 2 to 4), the solid fraction (6), the coordination number (7), and
  the mean stress `p` (9), close to the 1000 Pa imposed on the walls at the
  end. `WallBox 0 1 2 3 4 5`, in `state.txt`, gives the six walls.

  ```sh
  postpro state.txt
  gnuplot -p -e "plot 'SampleState.txt' using 1:6 with lines title 'solid fraction'"
  ```
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
   widens along `x` and `z`, and its solid fraction decreases (it dilates):
   `sh solidFraction.sh conf100` (run in the new folder, with a copy of the
   script) gives about 0.51 at `t = 1.5 s`.
   The tutorial [`triaxialTest`](../triaxialTest) does this test completely,
   with the stress-strain curves.

### Going further

- [`prepro_cpptools`](../prepro_cpptools): more about `generator`.
- [`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic): the same compression
  in a periodic cell, without walls.
- The other servos:
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html).
