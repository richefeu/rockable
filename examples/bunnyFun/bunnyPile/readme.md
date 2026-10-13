## 15. bunnyFun 1: a pile of bunnies

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

Sixty Stanford bunnies, of sizes between 3 and 4.7 cm, are stacked in a tube
by `generator`, without any overlap. They fall to the bottom of the tube and
settle, then the tube is lifted and the pile spreads on the floor. On the way,
the whole chain of pre-processing is used on a shape that is not convex: from
the mesh to the shape, from the shape to a packing, from the shape to the skin
drawn by `see`.

### Notions

- The converter **`stl2shape`**, and the name of the shape it writes
  ([converters](https://richefeu.github.io/rockable/converters.html)).
- **`shapeSurvey -P`**: the mass properties and the **oriented bounding box**
  (OBB) of a shape ([shapeSurvey](https://richefeu.github.io/rockable/shapeSurvey.html)).
- The packing command **`generatePacking:OBBs`**: any shape, through its
  bounding box, with a range of sizes
  ([generator](https://richefeu.github.io/rockable/generator.html)).
- A shape library made of several files, an **open surface** (the tube), and a
  body turned by its quaternion.
- **Stabilising** a packing: the numerical damping switched on and off with
  `Tempo` ([syntaxConf](https://richefeu.github.io/rockable/syntaxConf.html),
  [dissipation](https://richefeu.github.io/rockable/dissipation.html)).
- The skins of `see` for a whole library, with **`shape2mesh -s`**.

### Files

| File                | Content                                                              |
|---------------------|----------------------------------------------------------------------|
| `../bunny.stl`      | The Stanford bunny, 309 triangles, about 10 cm long (see [`../readme.md`](../readme.md)) |
| `genShapes.txt`     | Commands of `generator` that write `walls.shp`: the floor and the tube |
| `genPile.txt`       | Commands of `generator` that write `input.txt`: the bodies and the bunnies |
| `drivingSystem.txt` | The servo that lifts the tube                                        |
| `pile.sh`           | The number of bunnies, their heights, their distance to the axis, the solid fraction |

### Run it

```sh
stl2shape -r 0.001 -c ../bunny.stl               # writes mesh.shp
sed 's/^name .*/name Bunny/' mesh.shp > bunny.shp # names the shape Bunny
shapeSurvey bunny.shp -P bunny.shp               # mass properties and OBB
generator genShapes.txt                          # writes walls.shp
cat walls.shp bunny.shp > shapes.shp             # the library of shapes
generator genPile.txt                            # writes input.txt
rockable input.txt -j 4                          # about a minute, conf0 to conf90
shape2mesh shapes.shp --rmsh -s                  # the skins, shapes.rmsh
see
sh pile.sh conf0     # stacked by generator
sh pile.sh conf35    # settled, before the tube is lifted
sh pile.sh conf90    # spread on the floor
```

### What happens

**The shape.** `stl2shape` turns the mesh into a sphero-polyhedron of radius
1 mm (157 vertices, 464 edges once the 177 duplicated ones are removed by
`-c`, 309 faces). It names the shape after the file, here `../bunny.stl`,
which is not a convenient name: `sed` renames it `Bunny`. `shapeSurvey -P`
then computes the volume (136 cm³), the inertia, and the **OBB**, the smallest
box aligned with the principal axes that holds the whole shape, radius
included:

```
obb.extent 5.2513679405519001e-02 3.0797885114754682e-02 4.4848872436893027e-02
obb.center -4.3813175696563768e-03 4.0569623602885740e-03 4.0685712852794459e-03
```

The box is 10.5 × 6.2 × 9.0 cm (the extents are half-sizes), and its centre is
not the mass centre of the bunny: it is 7 mm away, since the bunny is not
symmetric. The bunny fills only 23% of its box.

**The packing.** In `genPile.txt`, the line

```
generatePacking:OBBs shapes.shp Bunny CYLINDER STACK Y  -0.0585 0 -0.0585  0.117 0.8 0.117  0.3 0.45  400 0.6  0 2  0
```

reads the OBB of `Bunny` in `shapes.shp`, and stacks boxes in a cylinder of
axis `Y` (the bounding box of the cylinder starts at `(-0.0585, 0, -0.0585)`
and is `0.117 × 0.8 × 0.117`, a little narrower than the tube). Each box gets
a homothety between 0.3 and 0.45 and a random orientation. In the mode
`STACK`, it is tried 20 times at random places, goes down to the lowest place
where it overlaps no other box, and the lowest trial is kept. Each box is then
replaced by its bunny, placed so that its box is where it was packed. The
packing stops at 400 bunnies, or when the boxes fill 60% of the cylinder, or,
as here, when no box finds a place:

```
@generatePacking:OBBs, no place found for box 62, the container is full
@generatePacking:OBBs, 61 particles, height 0.791285, solid fraction of the boxes = 0.234378, of the particles = 0.0550084
```

The boxes fill only 23% of the cylinder: they are placed one by one, and do
not move once placed. Since a bunny fills 23% of its box, the bunnies fill
only 5.5% of the cylinder. No bunny touches another one: the computation
starts without any overlap, but far from a dense state.

**The container.** `generateShape:thin_cylinder` writes a tube made of 24 flat
faces and no bottom (an **open surface**, `isSurface`), whose axis is `x`. In
`genPile.txt`, the tube is given the quaternion `0.7071 0 0 0.7071`, a quarter
turn about `z`, which brings its axis along `y`.

**The settling.** The bunnies fall to the bottom of the tube. A fall of 80 cm
would bring them at nearly 4 m/s, and the contacts, with a stiffness of
10⁴ N/m, would sink about 2 mm into the shapes, much more than their radius
(0.3 to 0.45 mm once scaled): a vertex could then cross a face. The line

```
Tempo NDCoeff Range 0 0.6 0.7 0
```

sets a strong numerical damping (0.7) until `t = 0.6 s`, and none after. It
slows the fall down, and the pile is nearly at rest at `t = 0.58 s`; when
the damping stops, the bunnies move a little, then settle again before the
tube is lifted (exercise 6 shows what happens without the damping). This is what **stabilising** a packing means: it
was made without overlap, but it is not in equilibrium, and a computation has
to bring it there.

**The collapse.** `drivingSystem.txt` lifts the tube at 0.5 m/s from
`t = 0.7 s` (`Servo ramp _y_Vel_ 1 0.0 0.5 0.7 0.71`), and the pile spreads on
the floor. `pile.sh` gives:

```
                         conf0      conf35     conf90
highest centre (m)       0.774      0.149      0.040
mean height (m)          0.382      0.066      0.014
mean distance to axis    0.028      0.036      0.083
solid fraction (tube)    0.053      0.278      -
```

The settling packs the bunnies five times more densely than `generator` did.
The pile is 15 cm high and 12 cm wide; once the tube is lifted, it spreads
into a layer of one or two bunnies, 4 cm high at most, and the bunnies end at
8.3 cm from the axis on average. The solid fraction in the tube is meaningless
once the tube is gone.

### What to look at

- In `see`, the tube hides the bunnies that are inside: the key `d` hides the
  floor and the tube. The keys `+` and `-` go through the conf-files.
- `see` finds `shapes.rmsh` next to `shapes.shp`, and draws each body with its
  skin: the bunnies have ears. Without `-s`, `shape2mesh` gives the tube a skin
  of a million triangles (91 MB); with `-s`, the flat parts are merged, and the
  file holds 145 000 triangles for the tube, 8 300 for the bunny (7.7 MB).
- The kinetic energy (`gnuplot checkplots.txt`): the fall, slowed by the
  damping, the rest, and the collapse. It does not fall to zero at the end: the
  tube still goes up.
- Other seeds give other packings: with `1` or `2` as the last value of the
  line `generatePacking:OBBs`, 61 and 59 bunnies, and a distance to the axis
  of 8.0 and 8.5 cm at the end.

### Exercises

1. **Spheres.** Add a ball of the same volume as the bunny to `genShapes.txt`,
   before `close`: `generateShape:sphere Ball 0.0319`. Run
   `generator genShapes.txt`, `cat walls.shp bunny.shp > shapes.shp`, write
   `Ball` instead of `Bunny` in the line `generatePacking:OBBs`, run
   `generator genPile.txt` and the computation (a few seconds). A ball fills
   52% of its box: 153 balls fit in the tube, and they settle at a solid
   fraction of 0.46. Once the tube is lifted, they roll away and fall off the
   floor: a ball has nothing that stops it rolling.
2. **Rolling resistance.** In the `input.txt` of the balls, add after each
   line `muContact`, for the pairs `0 0` and `0 1`:

   ```
   krContact 0 0 0.5
   murContact 0 0 0.005
   ```

   The balls now stop, in a heap 6.7 cm high, at 9.5 cm from the axis on
   average. A bunny needs no rolling resistance: its shape stops it.
3. **A mixture.** Write `Bunny,Ball` in the line `generatePacking:OBBs`: each
   grain is a bunny or a ball, at random. 87 grains fit in the tube, and the
   balls roll away when it is lifted.
4. **The mode `RANDOM`.** Write `RANDOM` instead of `STACK`: the boxes are put
   anywhere in the cylinder, the largest first. Only 42 bunnies fit: the boxes
   jam before they fill the cylinder. They settle more loosely (0.26).
5. **A wider range of sizes.** Write `0.25 0.5` instead of `0.3 0.45`: 55
   bunnies, from 2.6 to 5.3 cm, which settle a little more densely (0.30).
6. **Without damping.** Remove the line `print Tempo NDCoeff ...` of
   `genPile.txt`, set the friction between the bunnies to 0.3
   (`print muContact 0 0 0.3`), and run it again. The bunnies hit the pile
   at full speed; at `t = 0.34 s`, the kinetic energy of rotation jumps to
   3 J: one bunny has crossed the wall of the tube. It stays caught in the
   wall, and goes up with the tube, to 82 cm (`pile.sh conf90`, highest
   centre). Without the damping, the stiffness and the time step would have to
   be chosen for impacts at 4 m/s.

### Going further

- [`complexShape`](../../complexShape): the conversion of a scanned rock, and
  the options of `shapeSurvey`.
- [`columnCollapse`](../../columnCollapse): the collapse of a column of grains
  in a channel, and its run-out.
- All the options of `generatePacking:OBBs` (a cuboid container, a container
  without top, several shapes, the seed):
  [generator](https://richefeu.github.io/rockable/generator.html).
