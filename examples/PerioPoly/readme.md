## 8. A periodic sample of polyhedra

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

512 rounded dodecahedra, on a grid in a periodic cell, are compressed under a
pressure of 100 until the sample is dense. It is the same loading as
[`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic), with angular
particles: they pack differently, and they cost much more to compute.

### Notions

- **Sphero-polyhedra**: a shape made of vertices, edges and faces, dilated by
  a Minkowski radius ([syntaxShapes](https://richefeu.github.io/rockable/syntaxShapes.html)).
- The **cost** of the contact detection between polyhedra.
- The options of the periodic cell: `mh`, `cellVelocityCorrection`
  ([periodicBoundaryConditions](https://richefeu.github.io/rockable/periodicBoundaryConditions.html)).
- The **variability** of the final state of a granular sample.

### Files

| File                | Content                                                          |
|---------------------|------------------------------------------------------------------|
| `input.txt`         | The conf-file: the cell and 512 dodecahedra                       |
| `shape.shp`         | The shapes `POLYR0` (a rounded dodecahedron) and `SPHERE0`        |
| `drivingSystem.txt` | `PeriodicLoading IsotropicCompression 100`                        |
| `solidFraction.sh`  | Computes the size of the cell and the solid fraction              |
| `state.txt`         | Command file of `postpro`: the same, and more                     |

### Run it

```sh
rockable input.txt           # about 4 minutes
see conf30
sh solidFraction.sh
```

### What happens

**The shape.** `POLYR0` is a regular dodecahedron (20 vertices, 30 edges,
12 pentagonal faces) with a Minkowski radius of 0.055: its edges and corners
are rounded. Its size is about 1 (the units of this example have no
dimension: `density 0 0.2`, `knContact 0 0 1e5`, a pressure of 100). Its
volume and inertia are given in `shape.shp` (`preCompDone y`).

**The sample.** The 512 particles are on a grid of 8 × 8 × 8 in a cubic cell
of side 8, with random orientations: the solid fraction is 0.348 at the start.
There is no gravity, and the particles are all in group 0.

**The cell.** As in `isoComprSpheresPeriodic`, the cell is driven by a
pressure (`drivingSystem.txt`), its degrees of freedom have a mass set by
`mh .05`, and `numericalDampingCoeff 0.1` helps the sample reach its
equilibrium. `cellVelocityCorrection 2.5e-05`, equal to the time step,
removes at every step the mean velocity of the particles: with no wall, nothing
else would prevent the whole sample from drifting.

The cell shrinks from 8 to about 6.6–7.0 in 0.15 s, then hardly moves. The
solid fraction reaches about **0.58**.

**A variable result.** Run the computation twice with a different `interConf`:
the final sizes of the cell and the solid fraction are not exactly the same
(0.578 or 0.587, for instance). Each save of a conf-file converts the positions
from the coordinates of the cell to real coordinates and back, which changes
the last digits, and the dynamics of a granular sample amplifies such tiny
differences. This is not an error: a real sample would not give the same
value twice either. A result has to be checked on several samples.

### What to look at

- In `see`, the key `P` shows the cell. Compare `conf0`, `conf10` and
  `conf30`.
- `sh solidFraction.sh` prints the time, the three sizes of the cell and the
  solid fraction. It reads the volume of each shape in `shape.shp`, and
  multiplies it by the cube of the homothety of each particle.
- The same with `postpro`: `postpro state.txt` runs the post-processor
  `SampleState` (see [`postProcessing`](../postProcessing)) and writes
  `SampleState.txt`, with also the coordination number (column 7, about 6.7)
  and the mean stress (column 9, close to 100).
- `perf.txt` (second column): the number of seconds of simulation computed per
  second of real time. Compare it with that of `isoComprSpheresPeriodic`.

### Exercises

1. Replace the dodecahedra by spheres of the same volume: the shape `SPHERE0`,
   already in `shape.shp`, has a radius of 0.436.

   ```sh
   sed -i.bak 's/^POLYR0 /SPHERE0 /' input.txt
   ```

   The computation takes 2 seconds instead of 4 minutes, and the solid
   fraction reaches 0.61: the spheres pack more densely than the angular
   particles.
2. Run the computation with 4 threads (`rockable input.txt -j 4`): with 512
   polyhedra, it is about twice as fast.
3. Shear the compressed sample: copy `conf30` and `shape.shp` into a new
   folder, set `tmax 0.6` in `conf30`, write
   `PeriodicLoading SimpleShearDeformable XY 100 0.5` in `drivingSystem.txt`,
   and run `rockable conf30 -j 4` (a few minutes). The term `xy` of `h` (the
   second number of the line `h` in the conf-files) grows by 0.5 × `h_yy`,
   about 3.3, per second, while the three sizes are kept at a pressure of
   100.

### Going further

- [`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic): the same loading
  with spheres.
- The options of the periodic cell for large shears (`cellRelattice`,
  `cellDriveCauchy`):
  [periodicBoundaryConditions](https://richefeu.github.io/rockable/periodicBoundaryConditions.html).
