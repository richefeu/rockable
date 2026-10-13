## 3c. A triaxial test: compressing, restarting, shearing

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

The classical test of soil mechanics, in two computations. A sample of 343
spheres is first compressed isotropically at 1000 Pa, as in
[`isoComprSpheres`](../isoComprSpheres). The last conf-file of this
compression then becomes the initial state of a second computation, where
the top wall comes down at a constant velocity while the lateral walls keep
the pressure at 1000 Pa. The sample first resists more and more, reaches a
peak of strength, dilates, and ends in a steady state: the *critical state*.

### Notions

- **Restarting a computation** from a conf-file, and changing some of its
  parameters on the way
  ([commandLine](https://richefeu.github.io/rockable/commandLine.html)).
- The servo **`tritriBiaxialCompression`**: a velocity on one wall, a stress
  on the lateral ones
  ([servoControllers](https://richefeu.github.io/rockable/servoControllers.html)).
- The **mass of the walls** driven by a force.
- The **stress-strain curves** of a granular material, from `SampleState`
  ([postProcessing](https://richefeu.github.io/rockable/postProcessing.html)):
  peak and critical state, contraction and dilatancy.

### Files

There is one folder per computation, because the loading is always read in a
file named `drivingSystem.txt`, in the folder of the computation.

| File                          | Content                                                          |
|-------------------------------|------------------------------------------------------------------|
| `compression/input.txt`       | The conf-file of the compression: 6 walls and 343 spheres         |
| `compression/shapes.txt`      | The shapes `Sphere` and `x-wall`, `y-wall`, `z-wall`              |
| `compression/drivingSystem.txt` | The isotropic compression at 1000 Pa                            |
| `compression/genShapes.txt`, `compression/genSample.txt` | Commands of `generator`              |
| `shear/restart.sh`            | Builds `shear/input.txt` from the last conf-file of the compression |
| `shear/shapes.txt`            | The same shapes                                                   |
| `shear/drivingSystem.txt`     | The triaxial loading                                              |
| `shear/state.txt`             | The command file of `postpro`                                     |
| `shear/triaxial.sh`           | Strains and stresses, from `SampleState.txt`                      |

### Run it

```sh
cd compression
rockable input.txt           # about 8 seconds, conf0 to conf30
cd ../shear
sh restart.sh                # writes input.txt from ../compression/conf30
rockable input.txt           # about 30 seconds, conf0 to conf50
postpro state.txt            # writes SampleState.txt
sh triaxial.sh > triaxial.txt
```

### What happens

**The compression.** It is built as in `isoComprSpheres` (see its readme):
`generator` places the walls of a box of 7 cm and 7 × 7 × 7 spheres on a grid,
and `homothetyRange` gives the spheres a radius between 3 and 4.5 mm. Two
choices differ:

- There is **no friction** during the compression (`muContact 1 1 0.0`): the
  spheres slide on each other, and the sample gets dense (solid fraction
  0.567).
- The **walls are light** (`density 0 50`). A wall driven by a force is a body
  with a mass, which moves under the difference between the imposed force and
  the force of the grains. A wall heavier than the sample (with
  `density 0 2700` it weighs 0.4 kg, twice the sample) oscillates, and the
  stress it imposes oscillates with it (exercise 1).

**The restart.** A conf-file holds everything a computation needs: the
parameters, the bodies, and their contacts with the forces they carry. Any
conf-file can therefore be given to `rockable`, and the computation continues
from that state. `restart.sh` copies `conf30` of the compression into
`input.txt`, and changes a few lines with `sed`:

```
t 0                      the time starts again from 0
iconf 0                  so do the numbers of the conf-files
tmax 1.0                 the duration of the shear
interConf 0.02           a conf-file every 0.02 s, 51 in all
muContact 1 1 0.5        friction between the spheres
```

The files that are not in the conf-file, `shapes.txt`, `drivingSystem.txt` and
`dataExtractors.txt`, are read in the folder where `rockable` runs. It is the
`drivingSystem.txt` of the folder `shear` that changes the loading:

```
Servo tritriBiaxialCompression 0 1 2 3 4 5
1000 0.01
```

the top wall (`ymax`) comes down at 0.01 m/s, the walls `xmax` and `zmax` are
pushed at 1000 Pa, and the three `min` walls stay in place. The friction
between the spheres and the walls stays at zero: the top and bottom walls
must not hold the sample back as it widens (exercise 2).

**The results.** `triaxial.sh` gives, from `SampleState.txt`, the axial strain
`eps_y = 1 - Ly/Ly0`, the volumetric strain `eps_v = 1 - V/V0` (negative when
the sample dilates), the mean stress `p`, the deviatoric stress
`syy - (sxx + szz)/2`, the stress ratio `q/p`, the solid fraction and the
coordination number:

```
# time  eps_y  eps_v    p      syy-(sxx+szz)/2  q/p    solid_fraction  Z
0.000   0.0000  0.0000   998.9     1.6         0.0714  0.5666          5.07
0.200   0.0386 -0.0123  1454.0  1295.1         0.8994  0.5598          3.19
0.500   0.0965 -0.0344  1299.5  1137.0         0.8781  0.5478          3.32
1.000   0.1929 -0.0442  1367.0   942.2         0.6940  0.5426          3.40
```

- The **stress ratio** `q/p` jumps to about 0.75 at once (a dense sample is
  stiff), stays around 0.85 up to 10% of axial strain (the **peak**), then
  decreases to about 0.6: the **critical state**, where the sample deforms
  without changing its strength nor its volume.
- The **volume** increases from the start: the dense sample **dilates**
  (4.4% at the end), because the grains must climb over each other
  to move. The solid fraction goes from 0.567 to 0.543, and stays there after
  15% of axial strain.
- The **lateral stresses** `sxx` and `szz` (columns 11 and 13 of
  `SampleState.txt`) stay close to the 1000 Pa imposed by the servo, which
  checks that the loading is slow enough to be quasi-static.
- The curves are noisy: with 343 spheres, a few contacts that slip change the
  stress by 10%. A real study uses thousands of grains, or averages several
  samples.

With `q/p = M`, the friction angle of the material is given by
`sin φ = 3M / (6 + M)`: about 22° at the peak and 16° at the critical state,
which is low because spheres roll easily (exercise 4).

### What to look at

- In `see`, in the folder `shear`, the sample shortens and widens (the key `d`
  hides the walls).
- The stress ratio and the volumetric strain against the axial strain:

  ```sh
  gnuplot -p -e "plot 'triaxial.txt' using 2:6 with lines title 'q/p'"
  gnuplot -p -e "plot 'triaxial.txt' using 2:3 with lines title 'eps_v'"
  ```
- The lateral stresses, against the 1000 Pa of the servo:

  ```sh
  gnuplot -p -e "plot 'SampleState.txt' using 1:11 with lines title 'sxx', '' using 1:13 with lines title 'szz', 1000"
  ```

### Exercises

1. Make the walls heavy: `density 0 2700` in `shear/input.txt` (or in
   `restart.sh`, with one more `-e` line). The lateral stresses still average
   1000 Pa, but now swing between 300 and 1800 Pa, and so do the curves.
2. Give the walls a friction (`muContact 0 1 0.5`): the top and bottom walls
   hold the sample back as it widens, the lateral stresses rise to about
   1300 Pa although the servo still asks for 1000 Pa, and the sample looks
   stronger (`q/p` about 1.05). In the laboratory, the ends of the sample are
   lubricated for the same reason.
3. Prepare a **loose** sample: keep the friction during the compression
   (`muContact 1 1 0.5` in `compression/input.txt`), then run `restart.sh`
   and the shear again. The solid fraction starts at 0.539, the sample
   **contracts** first, `q/p` increases slowly, without a peak, and the sample
   ends at the same solid fraction, about 0.54, as the dense one: the critical
   state does not depend on the initial state.
4. Add a rolling resistance between the spheres, after the line
   `muContact 1 1 0.5` of `shear/input.txt`:

   ```
   krContact 1 1 4
   murContact 1 1 0.0004
   ```

   (in the law `Default`, `murContact` is a length, here a tenth of the
   radius). The sample is much stronger (`q/p` about 1.5 at the peak, 1.15 at
   the critical state) and dilates three times more.
5. **Resume** a computation and change its loading. Copy the folder `shear`
   after the computation, and in the copy, write `1000 -0.01` in
   `drivingSystem.txt`, then run `rockable conf25`. The computation continues
   from `t = 0.5` (and writes `conf26`, `conf27`... over the old ones): the
   top wall goes up, the deviatoric stress falls to zero at once (the sample
   unloads elastically), then becomes negative: the sample is now in
   extension. It contracts first, then dilates again when `q/p` reaches about
   0.7 in extension.
6. `drivingSystem.txt` is read again at each conf-file: edit it **during** the
   computation, for instance to change the confining pressure, and see the
   effect in the next conf-files.

### Going further

- [`postProcessing`](../postProcessing): `SampleState`, the forces and
  ParaView.
- [`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic) and
  [`PerioPoly`](../PerioPoly): samples without walls, in a periodic cell, which
  avoids the effect of the walls.
- The other servos (`tritriCustom` for any combination of velocities and
  stresses):
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html).
