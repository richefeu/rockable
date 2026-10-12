## 7. Periodic boundary conditions: isotropic compression

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

The same 125 spheres as in [`isoComprSpheres`](../isoComprSpheres), compressed
the same way, but in a **periodic cell** instead of a box of walls. A sphere
that leaves the cell through one face comes back through the opposite face,
and it touches the images of the spheres on the other side: there is no wall,
and no wall effect.

### Notions

- The **periodic cell**: `usePeriodicCell`, the cell matrix `h`, its velocity
  `vh` and acceleration `ah`, the mass of its degrees of freedom `mh`
  ([periodicBoundaryConditions](https://richefeu.github.io/rockable/periodicBoundaryConditions.html)).
- The loading of the cell, **`PeriodicLoading`**, written in
  `drivingSystem.txt` in place of the `Control` and `Servo` entries used for
  walls.
- The **numerical damping**, `numericalDampingCoeff`.
- Reading the conf-files with a script (`h` and the homotheties).

### Files

| File                | Content                                                      |
|---------------------|--------------------------------------------------------------|
| `input.txt`         | The conf-file: the cell and 125 spheres                       |
| `shapes.txt`        | The shape `Sphere` (unit radius) and the shapes of walls (unused) |
| `drivingSystem.txt` | The loading of the cell                                       |
| `solidFraction.sh`  | Computes the size of the cell and the solid fraction          |

### Run it

```sh
rockable input.txt           # about 6 seconds
see conf30
sh solidFraction.sh          # the solid fraction during the compression
```

### What happens

**The cell.** In `input.txt`,

```
usePeriodicCell 1
h   0.05 0 0    0 0.05 0    0 0 0.05
vh   0 0 0    0 0 0    0 0 0
ah   0 0 0    0 0 0    0 0 0
mh 1
```

activates the periodic cell and makes it a cube of 5 cm. The columns of the
matrix `h` are the three vectors that span the cell (a sheared cell has
off-diagonal terms). The cell is not a fixed frame: its nine components are
degrees of freedom, integrated along with the particles, with a mass set by
`mh`. This is what allows a stress to be imposed on it. There is no wall in the
list of particles, and `nDriven 0`. The parameters of the group pair `0 1`
(`# sphere-wall`) are not used here: all the spheres are in group 1.

**The loading.** `drivingSystem.txt` holds

```
PeriodicLoading IsotropicCompression 100
```

The three diagonal components of `h` are driven by a stress of -100 Pa (a
compression of 100 Pa), and the shear components are kept at zero. It is the
periodic counterpart of the servo `tritriIsostaticCompression` of
`isoComprSpheres`.

**The damping.** `numericalDampingCoeff 0.1` damps the accelerations of the
particles, so that the sample reaches its equilibrium faster. It is a
numerical tool, not a physical dissipation: it must be small, or switched off,
when the dynamics matters.

**`mh`.** If the cell is too light, it reacts too fast to the contact forces
and the computation diverges: with `mh 0.02`, the cell matrix is `nan`
at `t = 0.01 s`. With `mh 1`, the cell shrinks smoothly.

The cell goes from 50 mm to about 30 × 29 × 30 mm in 0.15 s, then hardly
moves. The solid fraction goes from 0.13 to about **0.61**. With walls,
`isoComprSpheres` only reached 0.52 (at a pressure ten times higher): with so
few spheres, the layers along the walls, where the spheres are badly packed,
are a large part of that sample. The periodic cell has no such layer.

### What to look at

- In `see`, the key `P` shows or hides the cell. The spheres are drawn at their
  real position: a sphere that straddles a face is drawn on one side only, and
  its image on the other side is not drawn.
- `sh solidFraction.sh` reads `t` and `h` in each conf-file, and the radius of
  each sphere (its homothety, since the shape has a unit radius). It prints
  the time, the three sizes of the cell and the solid fraction:

  ```sh
  sh solidFraction.sh > phi.txt
  gnuplot -p -e "plot 'phi.txt' using 1:5 with lines title 'solid fraction'"
  ```
- The line `h` of the conf-files: the off-diagonal terms stay at zero.

### Exercises

1. Set the friction between the spheres to zero (`muContact 1 1 0.0`): the
   sample is denser (solid fraction about 0.64).
2. Compress at 1000 Pa: the compression is faster, but the final solid
   fraction is almost the same (about 0.61). The spheres are stiff (`knContact
   1e8`): the packing is set by the friction, not by the pressure.
3. Set `mh 0.02` and watch the line `h` of the conf-files.
4. Use the compressed sample for a triaxial test. Copy `conf30`, `shapes.txt`
   and `solidFraction.sh` in a new folder, set `tmax 0.6` in `conf30`, and
   write this `drivingSystem.txt`:

   ```
   PeriodicLoading TriaxialCompression Y 100 -0.5
   ```

   then run `rockable conf30`. The size along `y` decreases at a strain rate of
   0.5 /s (a negative rate shortens the cell), while the two other directions
   stay at 100 Pa. `sh solidFraction.sh 30 60` shows that the sample first
   gets denser (0.612 around `t = 0.42 s`), then dilates (0.608 at `t = 0.6 s`).

### Going further

- [`PerioPoly`](../PerioPoly): a periodic sample of polyhedra.
- The other loadings (`TriaxialCompression`, `SimpleShearDeformable`) and the
  options of the cell:
  [periodicBoundaryConditions](https://richefeu.github.io/rockable/periodicBoundaryConditions.html).
