## 14. Collapse of a granular column

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

540 grains fall into a narrow reservoir closed by a gate, and form a column
twice as high as it is long. At `t = 0.6 s` the gate is lifted: the column
collapses and spreads along a channel, and stops as a deposit with a gentle
slope. How far it goes depends on the shape of the column, on the friction,
and on the force law: with the law `GeoVisc`, it never stops.

### Notions

- A **time-dependent loading** with the servo **`ramp`**: a body that waits,
  then moves ([servoControllers](https://richefeu.github.io/rockable/servoControllers.html)).
- The preparation and the test in **one computation**.
- A grain made of **four spheres** (a shape with 4 vertices and no edge).
- The **run-out** of a granular collapse, and the profile of the deposit,
  read from a conf-file with `awk`.
- The force law **`GeoVisc`**, whose friction is viscous
  ([forceLaws](https://richefeu.github.io/rockable/forceLaws.html)).

### Files

| File                | Content                                                          |
|---------------------|------------------------------------------------------------------|
| `input.txt`         | The conf-file: the floor, the walls, the gate and 540 grains     |
| `shapes.txt`        | The shapes `Floor`, `EndWall`, `SideWall` and `Grain`             |
| `genShapes.txt`     | Commands of `generator` that write the shapes of the walls       |
| `genColumn.txt`     | Commands of `generator` that write the lines of the grains       |
| `drivingSystem.txt` | The servo that lifts the gate                                    |
| `profile.sh`        | The height of the deposit along the channel, and its front       |

### Run it

```sh
rockable input.txt       # about 70 seconds, conf0 to conf65
see
sh profile.sh conf29 > column.txt     # the column before the gate is lifted
sh profile.sh conf65 > deposit.txt    # the deposit
```

### What happens

**The set-up.** The channel is 5 cm wide and 40 cm long. Its five bodies are
driven (`nDriven 5`): the floor, the two side walls, the end wall at `x = 0`,
and the gate at `x = 6 cm` (body 4). They are cuboids made by `genShapes.txt`.
The grains are all the same shape, `Grain`: four spheres of radius 3.1 mm
whose centres are about 3 mm apart, a rounded grain of about 1 cm with a few
bumps. It is written by hand in `shapes.txt` (4 vertices, no edge, no face).
`genColumn.txt` places 5 × 27 × 4 of them on a grid, with random orientations
(the last `1`), 1.1 cm apart vertically.

**The loading.** `drivingSystem.txt` holds a single line,

```
Servo ramp _y_Vel_ 4 0.0 1.0 0.6 0.61
```

The vertical velocity of the gate is 0 until `t = 0.6 s`, grows to 1 m/s
between 0.6 and 0.61 s, then stays at 1 m/s. The floor and the other walls
have no control: they do not move. So, in a single computation, the grains
first fall and settle behind the gate (0.6 s is enough), then the column is
released.

**The collapse.** `profile.sh` cuts the channel into slices of 1 cm, and gives
the height of the highest grain in each slice. The front of the deposit is
where the first empty slice starts: the few grains that roll ahead are left
out.

```
                         height H     length L    aspect ratio a = H0/L0
column (conf29)          11.6 cm      6 cm        1.93
deposit (conf65)          7.0 cm     23 to 24 cm
```

The run-out, `ΔL = L - L0`, is about 17.5 cm, that is `ΔL/L0 ≈ 2.9`. The top
of the column comes down to 61% of its height, and the deposit ends with a
thin layer of one or two grains. The flow lasts about 0.6 s: the front no
longer moves after `t = 1.2 s` (conf60).

Experiments on granular columns, in channels (Lajeunesse et al., 2005; Lube
et al., 2005), find that the relative run-out `ΔL/L0` depends mainly on the
aspect ratio `a`, and grows in proportion to it when the column is not much
higher than long. Exercise 1 checks it.

### What to look at

- In `see`, the key `d` hides the floor, the walls and the gate, which are in
  front of the grains. The keys `+` and `-` go through the conf-files.
- The column and the deposit:

  ```sh
  gnuplot -p -e "set size ratio -1; plot 'column.txt' with steps title 'column', 'deposit.txt' with steps title 'deposit'"
  ```
- The kinetic energy (`gnuplot checkplots.txt`): the fall of the grains, their
  settling, the collapse. It does not fall to zero at the end: the gate still
  goes up at 1 m/s.

### Exercises

1. A lower column: in `genColumn.txt`, write `0.06 0.143 0.05` and `5 13 4`
   (13 layers), run `generator genColumn.txt`, and put its 260 lines in place
   of the grains of `input.txt` (with `Particles 265`). The column is now as
   high as long (`H0 = 6.1 cm`, `a = 1.0`), and spreads to 15 cm:
   `ΔL/L0 = 1.5`, half the run-out of the column twice as high. Its top comes
   down to 88% of its height only.
2. Less friction between the grains (`muContact 1 1 0.3`): the deposit
   reaches 27 cm.
3. Lift the gate slowly, at 0.1 m/s (`Servo ramp _y_Vel_ 4 0.0 0.1 0.6 0.61`,
   and `tmax 2.0`): the gate releases the column gradually and holds back a
   part of it, so the deposit is higher (8.1 cm) and shorter (16 cm). In
   the experiments, the gate is lifted as fast as possible.
4. Replace the law by `GeoVisc`: `forceLaw GeoVisc`, and replace the two lines
   `ktContact` by `viscTContact 0 1 0.01` and `viscTContact 1 1 0.01` (a
   tangential viscosity, in N·s/m). Run it with `tmax 2.6`. The tangential
   force is now proportional to the sliding velocity: a grain at rest carries
   no friction. The column collapses as before, but the deposit never stops:
   it creeps until the grains fall off the end of the floor. This law is meant
   for slow landslides, where the friction is viscous, not for a deposit at
   rest.
5. Add a rolling resistance between the grains, after `muContact 1 1 0.5`:

   ```
   krContact 1 1 0.5
   murContact 1 1 0.0005
   ```

   (in the law `Default`, `murContact` is a length, here a tenth of the radius
   of a grain). The grains roll less: the column settles looser (13.1 cm
   high instead of 11.6), and the deposit is steeper, with a top at 9.0 cm,
   69% of the column.

### Going further

- [`rockfall`](../rockfall): blocks rolling and sliding down a slope, and the
  law `Avalanche`.
- The other servos (`shaker`, `sawtooth_shaker`...):
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html).
- Larger columns, in parallel: `rockable input.txt -j 4`
  ([commandLine](https://richefeu.github.io/rockable/commandLine.html)).
