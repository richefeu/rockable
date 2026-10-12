## 11. From an STL mesh to a shape

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

A rock scanned as a triangulated mesh (`rock.stl`) is converted into a
Rockable shape with `stl2shape`, checked and pre-computed with `shapeSurvey`,
then used in a computation: the rock spins freely in space.

### Notions

- The converter **`stl2shape`** and its Minkowski radius
  ([converters](https://richefeu.github.io/rockable/converters.html)).
- The **pre-computation** of a shape (`preCompDone`, volume, inertia), and the
  application **`shapeSurvey`**
  ([shapeSurvey](https://richefeu.github.io/rockable/shapeSurvey.html)).
- The **homothety**, to change the unit of a shape.
- A library of shapes made of several files put together.

### Files

| File          | Content                                                                 |
|---------------|-------------------------------------------------------------------------|
| `rock.stl`    | A scanned rock: a binary STL mesh of 24 236 triangles, in millimetres    |
| `cubo5mm.stl` | A cubic frame (a non-convex shape), in millimetres                       |
| `input.txt`   | The conf-file: the rock, spinning                                        |

### Run it

```sh
stl2shape -r 0.5 -c rock.stl     # writes mesh.shp
mv mesh.shp rock.shp
shapeSurvey rock.shp             # optional: look at it, press C then s
rockable input.txt               # about 20 seconds, less than 1 s once pre-computed
see
```

### What happens

**The conversion.** `stl2shape` takes the mesh as the skeleton of a
sphero-polyhedron: its vertices, edges and triangles become those of the
shape, and the whole is dilated by the Minkowski radius given with `-r` (here
0.5, in the unit of the mesh: 0.5 mm). `-c` removes the duplicated edges: an
STL file stores each triangle on its own, so that every edge appears twice
(14 867 of the 51 214 edges are removed here). The shape is written in
`mesh.shp`, and is named after the STL file:

```
<
name rock.stl
radius 0.5
preCompDone n
...
nv 12113
```

`stl2shape` always writes `mesh.shp`, hence the `mv`, before converting
another mesh.

**The pre-computation.** `preCompDone n` means that the volume, the position
of the mass centre, the inertia and the principal axes of the shape are not
known yet. Rockable computes them when the computation starts (by a
Monte-Carlo integration), and prints them:

```
Computation of mass properties (volume, mass-center, inertia and body-frame) for shape rock.stl
                                        Volume: 42878.3
                                  inertia/mass: 179.796 215.449 210.393
```

For this rock, it takes about 20 seconds, at each start. `shapeSurvey rock.shp`
shows the shape; in its window, the key `C` computes the mass properties, and
`s` saves them in `rock.shp` (which is overwritten), with `preCompDone y`: the
next computations start at once.

**The units.** The mesh is in millimetres: the volume is 42 878 mm³, about
43 cm³. In `input.txt`, the particle has a homothety of `1.0e-3`, which turns
the millimetres into metres; the radius, the mass and the inertia follow. With
a density of 2700 kg/m³, the rock weighs about 116 g.

**The computation.** There is no gravity and no other body. The rock is given
an angular velocity `(1, 6, 0.5)` rad/s, and spins freely for one second. Its
angular momentum is kept, but not its angular velocity: the components of
`vrot` in the conf-files vary along the rotation, because the inertia of the
rock is not the same about all its axes.

### What to look at

- In `see`, the rock and its rotation. The drawing of a mesh of this size is
  slow; `shape2mesh rock.shp --rmsh` writes a skin mesh, `rock.rmsh`, that `see`
  uses for a smoother drawing.
- The line of the rock in the conf-files: `Q` (the orientation) and `vrot`.
- The cost of a shape grows with its number of vertices, edges and faces. A
  mesh of 12 000 vertices is fine for one particle, much too heavy for a
  sample of thousands: decimate the mesh before converting it.

### Exercises

1. Convert the rock with a larger radius, `stl2shape -r 2 -c rock.stl`, and
   compare the two shapes in `shapeSurvey`: the details are smoothed, and the
   volume is larger.
2. Convert `cubo5mm.stl`, a cubic frame of 120 mm made of bars of 5 mm. The
   shape is not convex: a grain could go through it.
3. Let the rock fall on a plane. Put the shape `Plan` of
   [`helloworld`](../helloworld) and the rock in a single library:

   ```sh
   awk '/^</{b=$0; k=0; next} /^name Plan/{k=1; print b} k{print} /^>/{k=0}' ../helloworld/shapes.shp > plan.shp
   cat rock.shp plan.shp > shapes.shp
   ```

   then, in a copy of `input.txt`: `shapeFile shapes.shp`, `gravity 0 -9.81 0`,
   `dt 1e-5`, `interVerlet 0.001`, `nDriven 1`, `density 1 2700`, the contact
   parameters of the pair of groups `0 1` (`knContact 0 1 1e4`,
   `en2Contact 0 1 0.1`, `ktContact 0 1 1e4`, `muContact 0 1 0.5`, and `0` for
   `krContact` and `murContact`), and the two particles:

   ```
   Particles 2
   Plan     0 0 0.1      0 0 0     0 0 0   0 0 0   1 0 0 0   0 0 0   0 0 0
   rock.stl 1 0 1.0e-3   0 0.1 0   0 0 0   0 0 0   1 0 0 0   1.0 6.0 0.5   0 0 0
   ```

   The plane is scaled down to 40 × 10 cm. The rock lands at `t = 0.1 s`, and
   is at rest from about `t = 0.6 s`.

### Going further

- [`helloworldPlus`](../helloworldPlus): a converted shape in a small system.
- The options of `stl2shape` (`-s`, `-m`, `-z`) and the other converters:
  [converters](https://richefeu.github.io/rockable/converters.html).
- The format of the shapes:
  [syntaxShapes](https://richefeu.github.io/rockable/syntaxShapes.html).
