## 3b. Post-processing: `postpro`, the conf-files, and ParaView

### Goal

The sample of [`isoComprSpheres`](../isoComprSpheres) is compressed again,
then analysed after the computation, from its conf-files: solid fraction,
coordination number and stress with `postpro`, distribution of the contact
forces with a short script, and a 3D view of the particles and of the forces in
ParaView.

### Notions

- **`postpro`**, its command file, and the post-processor **`SampleState`**
  ([postProcessing](https://richefeu.github.io/rockable/postProcessing.html)).
- The **stress tensor** of a granular sample, computed from the contact forces.
- The block **`Interactions`** of the conf-files, read with `awk`
  ([syntaxConf](https://richefeu.github.io/rockable/syntaxConf.html)).
- **`conftovtk`** and **ParaView**
  ([visualisation](https://richefeu.github.io/rockable/visualisation.html)).

### Files

| File                   | Content                                                            |
|------------------------|--------------------------------------------------------------------|
| `input.txt`            | The conf-file of `isoComprSpheres`, with `tmax 0.3` and `interConf 0.005` |
| `shapes.txt`           | The shapes `Sphere` and the walls                                   |
| `drivingSystem.txt`    | The servo: a pressure of 1000 Pa on the walls                       |
| `state.txt`            | The command file of `postpro`                                       |
| `forceDistribution.sh` | The distribution of the normal contact forces in a conf-file       |

### Run it

```sh
rockable input.txt             # about 8 seconds, 61 conf-files
postpro state.txt              # writes SampleState.txt
sh forceDistribution.sh conf60 > fn.txt
conftovtk                      # writes the VTK files, for ParaView
```

### What happens

**The computation.** It is the isotropic compression of `isoComprSpheres`:
125 spheres in a box of six walls, the three `max` walls pushed by a pressure
of 1000 Pa. A conf-file is saved every 0.005 s, so that the analysis has 61
states to work on.

**`postpro`.** It reads the conf-files one after the other, and runs a
post-processor on each of them. Its command file, `state.txt`, is

```
firstConf 0
lastConf 60
stepConf 1

PostProcessor SampleState
WallBox 0 1 2 3 4 5
```

`SampleState` writes one line per conf-file in `SampleState.txt`. `WallBox`
gives the six walls of the box (`xmin xmax ymin ymax zmin zmax`): the sizes of
the sample are the distances between their inner faces. With a periodic cell,
this line is not needed: the cell gives the volume.

| Column | Content                                                           |
|--------|-------------------------------------------------------------------|
| 1      | time                                                              |
| 2-4    | sizes of the sample                                               |
| 5      | volume                                                            |
| 6      | solid fraction (volume of the particles over the volume)          |
| 7      | coordination number `Z` (mean number of contacts per particle)    |
| 8      | number of contacts                                                |
| 9, 10  | mean stress `p` and von Mises stress `q`                           |
| 11-16  | `sxx syy szz sxy sxz syz`                                          |

**The stress.** It is computed from the contact forces only, as
`(1/V) sum over the particles and their contacts of F (x) (x_p - x_c)`: the
force received at a contact times the vector from the contact point to the
centre of the particle. It is positive in compression. At the end of the
compression, `sxx`, `syy` and `szz` are close to the 1000 Pa imposed on the
walls, which checks both the servo and the computation:

```
t     solid fraction   Z      p
0     0.128            0      0
0.1   0.515            1.01   383
0.15  0.517            2.52   930
0.3   0.518            2.62   871
```

The stress fluctuates around 1000 Pa (between 870 and 1040 Pa): with 125
spheres, a few contacts more or less change it. The coordination number,
about 2.6, is low: many spheres carry no force (the *rattlers*), and the
friction allows loose packings.

**The contact forces.** `forceDistribution.sh` reads the block `Interactions`
of a conf-file (one line per interaction, the normal force in the 16th column),
and prints the probability density of `fn/<fn>`, the normal force divided by
its mean. Most contacts carry less than the mean force, and a few carry two or
three times more: the forces are carried by *force chains*.

**ParaView.** `conftovtk` converts every conf-file of the folder into VTK
files: `Shepres0000.vtk`, `Shepres0001.vtk`... for the spheres (with their
radius and velocity), `Forces####.vtk` for the contacts (at the contact
points, with the normal force `Fn` and the tangential force `Ft`), and
`Polyr####.vtk` for the walls. In ParaView:

1. *File > Open*, choose the group `Shepres..vtk`, then *Apply*.
2. *Filters > Glyph*, with *Glyph Type* `Sphere`, *Orientation Array* `No
   orientation array`, *Scale Array* `Radius`, *Scale Factor* `2` (the sphere of
   a glyph has a radius of 0.5), *Glyph Mode* `All Points`; *Apply*, and colour
   by `Vel`.
3. Open the group `Forces..vtk` the same way, and add a `Glyph` with *Scale
   Array* `Fn` and a *Scale Factor* of about `0.01`; colour it by `Fn`. Lower
   the opacity of the spheres to see the forces inside the sample.
4. The play button of ParaView goes through the conf-files.

### What to look at

- The solid fraction, the coordination number and the mean stress against
  time:

  ```sh
  gnuplot -p -e "plot 'SampleState.txt' using 1:6 with lines title 'solid fraction'"
  gnuplot -p -e "plot 'SampleState.txt' using 1:7 with lines title 'Z'"
  gnuplot -p -e "plot 'SampleState.txt' using 1:9 with lines title 'p (Pa)'"
  ```
- The distribution of the forces:

  ```sh
  gnuplot -p -e "set logscale y; plot 'fn.txt' using 1:2 with linespoints title 'pdf of fn/<fn>'"
  ```
- In ParaView, the largest forces (in red) and the chains they form.

### Exercises

1. Run `postpro` with a command file that processes one conf-file out of
   ten (`stepConf 10`).
2. Compare the solid fraction of `SampleState.txt` with the one of the script
   `solidFraction.sh` of `isoComprSpheres`: they are the same.
3. Use `SampleState` on [`isoComprSpheresPeriodic`](../isoComprSpheresPeriodic):
   write a command file without `WallBox` (the periodic cell gives the
   volume), and check that `p` is close to the 100 Pa imposed on the cell.
4. Process the conf-files of [`PerioPoly`](../PerioPoly): the dodecahedra
   reach a coordination number of about 6.7, much more than the spheres.
5. Run `forceDistribution.sh` on the last conf-file of `PerioPoly` (512
   particles, more contacts): the distribution is smoother.
6. The post-processor `ParticleStress` writes the tensorial moment of each
   particle in `ParticleStress-conf60.txt`, and prints their sum divided by
   the volume given with `Volume` (take it in the 5th column of
   `SampleState.txt`): it is the stress tensor of `SampleState`.

### Going further

- The other post-processors, and how to write a new one in C++:
  [postProcessing](https://richefeu.github.io/rockable/postProcessing.html),
  [developerGuide](https://richefeu.github.io/rockable/developerGuide.html).
- The data extractors, which measure during the computation instead of
  afterwards: [dataExtractors](https://richefeu.github.io/rockable/dataExtractors.html).
