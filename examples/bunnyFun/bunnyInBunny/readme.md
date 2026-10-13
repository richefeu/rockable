## 16. bunnyFun 2: bunnies in a bunny

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

A big hollow bunny, 26 cm long, has a hole in its back. Thirty small bunnies
fall into it through the hole. Then the big bunny rolls over, slowly, and
empties. The big bunny is made from the same mesh as the small ones: the shape
file written by `stl2shape` is edited with a small `awk` script, which cuts
the hole and turns the solid into a shell.

### Notions

- **Two shapes from one mesh**, with the options `-n` and `-o` of `stl2shape`
  ([converters](https://richefeu.github.io/rockable/converters.html)).
- **Editing a shape file**: removing faces, edges and vertices, and the
  keyword **`isSurface`**, which makes a shape a shell
  ([syntaxShapes](https://richefeu.github.io/rockable/syntaxShapes.html)).
- A **container that is not convex**, and the contacts on both sides of its
  faces.
- A **driven rotation**: the servo `ramp` on the control `_xrot_Vel_`
  ([drivingSystem](https://richefeu.github.io/rockable/drivingSystem.html),
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html)).
- The **orientation of a body** in a conf-file (its quaternion), used by a
  script that counts the bunnies inside the big one.

### Files

| File                | Content                                                                    |
|---------------------|----------------------------------------------------------------------------|
| `../bunny.stl`      | The Stanford bunny, 309 triangles (see [`../readme.md`](../readme.md))      |
| `hollow.sh`         | Turns a shape into a shell with a hole (sh + awk)                           |
| `genShapes.txt`     | Commands of `generator` that write `floor.shp`                              |
| `genInput.txt`      | Commands of `generator` that write `input.txt`: the floor, the big bunny, the small ones |
| `drivingSystem.txt` | The servo that turns the big bunny                                         |
| `count.sh`          | The number of small bunnies inside and outside the big one (sh + awk)       |

### Run it

```sh
stl2shape -r 0.001 -c -n Bunny -o bunny.shp ../bunny.stl      # the small bunny
shapeSurvey bunny.shp -P bunny.shp
stl2shape -r 0.001 -c -n BigBunny -o raw.shp ../bunny.stl     # the big one
sh hollow.sh raw.shp big.shp 0.006 0.026 -0.01 0.021          # a shell with a hole
generator genShapes.txt                                       # writes floor.shp
cat floor.shp big.shp bunny.shp > shapes.shp
generator genInput.txt                                        # writes input.txt
rockable input.txt -j 4                                       # 1.5 minutes, conf0 to conf95
shape2mesh shapes.shp --rmsh -s
see
for i in 0 35 70 75 80 85 90 95; do sh count.sh conf$i; done
```

### What happens

**The shell.** The mesh is converted twice. `Bunny` is the small bunny, a
solid, pre-computed by `shapeSurvey` as in [`bunnyPile`](../bunnyPile).
`BigBunny` is written in `raw.shp`, and `hollow.sh` turns it into a shell:

```
sh hollow.sh raw.shp big.shp 0.006 0.026 -0.01 0.021
33 faces removed, 276 kept; 147 vertices, 422 edges
```

The script reads the shape file (its format is plain text: the vertices after
`nv`, the edges after `ne`, the faces after `nf`), removes the 33 faces whose
centre is closer than 2.1 cm to the point `(0.006, 0.026, -0.01)`, on the
back, then the edges and the vertices that no longer belong to any face (they
would stay across the hole, as rods and balls). It adds the keyword
`isSurface`: the faces now make a shell of thickness twice the radius, and a
small bunny can touch them from both sides. Last, it writes `preCompDone y`,
with a box that bounds the shape, the volume of the shell and a rough
inertia: the big bunny is a driven body, whose mass properties hardly matter.

The shape is not given to `shapeSurvey`: the coordinates of the hole are
those of the mesh, and `shapeSurvey` would move the shape into its principal
frame. In `genInput.txt`, the big bunny gets a homothety of 2.5 (26 cm long,
a hole of about 7 cm), and its centre is 20 cm above the floor, so that it
can turn without touching it.

**The filling.** `generatePacking:OBBs` stacks 31 small bunnies (2.6 to
3.7 cm) in a cylinder of 5.5 cm diameter above the hole. As in `bunnyPile`,
a strong numerical damping (`Tempo NDCoeff Range 0 1.4 0.7 0`) slows their
fall. They go through the hole and settle in the big bunny; two of them
bounce on the edge of the hole and fall on the floor.

**The emptying.** `drivingSystem.txt` holds

```
Servo ramp _xrot_Vel_ 1 0.0 -1.5 1.6 1.7
```

the angular velocity of the big bunny (body 1) about `x`, its axis from the
nose to the tail, goes from 0 to -1.5 rad/s between 1.6 and 1.7 s, then stays
there. A driven body turns about its own centre (its position in the
conf-file, the origin of its shape). The big bunny rolls on its side, its
back turns down, and the small bunnies pour out through the hole.

**Counting.** `count.sh` reads the shape `BigBunny` in `shapes.shp`, and, in a
conf-file, the position, the homothety and the quaternion of the big bunny.
It brings the centre of each small bunny into the frame of the big one, and
casts a ray from it along `-y` (of that frame): the ray crosses the shell an
odd number of times when the small bunny is inside. The hole is on the back,
above the inside, so that the rays that go down from the inside do not go
through it.

```
t = 0.00 s, big bunny turned by   0 degrees: 0 small bunnies inside, 31 outside
t = 1.40 s, big bunny turned by   0 degrees: 29 small bunnies inside, 2 outside
t = 2.80 s, big bunny turned by  99 degrees: 29 small bunnies inside, 2 outside
t = 3.00 s, big bunny turned by 116 degrees: 28 small bunnies inside, 3 outside
t = 3.20 s, big bunny turned by 133 degrees: 24 small bunnies inside, 7 outside
t = 3.40 s, big bunny turned by 150 degrees: 13 small bunnies inside, 18 outside
t = 3.60 s, big bunny turned by 168 degrees: 6 small bunnies inside, 25 outside
t = 3.80 s, big bunny turned by 175 degrees: 5 small bunnies inside, 26 outside
```

Nothing comes out before the big bunny has turned by about 110 degrees: the
hole is on the back, and the small bunnies stay in the belly until the back is
well below the horizontal. Then most of them pour out within 0.4 s, between 130
and 170 degrees. Five of them are still inside at the end, in a recess of
the belly, on the side of the head, next to the hole.

### What to look at

- In `see`, the shell of the big bunny is transparent: the small bunnies are
  seen inside (the keys `a` and `A` change the transparency).
- `shapes.rmsh`, the skins: 139 000 triangles for the shell of the big bunny,
  which `shape2mesh` meshes on both sides.
- In the conf-files, the line of `BigBunny`: its quaternion `Q` changes, its
  position does not.

### Exercises

1. **A smaller hole.** Run `hollow.sh` with a radius of `0.016` instead of
   `0.021`, then `cat`, `generator genInput.txt` and the computation again (the
   column of small bunnies does not change). The hole is about 5 cm wide: 28
   bunnies get in, three bounce on its edge. When the big bunny turns, they
   come out one by one, and 19 are still inside at 175 degrees.
2. **Turn about another axis.** In `drivingSystem.txt`, write `_zrot_Vel_`
   instead of `_xrot_Vel_`: the big bunny now turns about `z`, its axis from
   side to side, and dips its head. The hole goes down much later: at 175
   degrees, 24 bunnies are still inside. Make `tmax` longer to see them come
   out.
3. **Less friction** between the small bunnies (`print muContact 0 0 0.3` in
   `genInput.txt`): they slide out a little earlier (21 inside at 133 degrees
   instead of 24), and only 4 stay inside.
4. **Without damping.** Remove the line `print Tempo NDCoeff ...`. The small
   bunnies fall faster, and here nothing goes wrong: all 31 get in, and the
   kinetic energy of rotation stays small. In [`bunnyPile`](../bunnyPile), the
   same change let a bunny cross a thin wall: the damping is a precaution,
   whose need depends on the speeds, the stiffness and the sizes.
5. **Move the hole.** The 4th to 6th arguments of `hollow.sh` are the centre
   of the hole, in the frame of the mesh (in metres, the bunny is about 10 cm
   long; `y` is up, the nose is towards `-x`). Open the hole elsewhere, and
   place the column of small bunnies above it: the column is a cylinder whose
   bounding box starts at `(x - 0.0275, 0.30, z - 0.0275)`, where `(x, z)` is
   the centre of the hole multiplied by 2.5.

### Going further

- [`bunnyPile`](../bunnyPile): the packing with `generatePacking:OBBs`, and the
  stabilisation of a packing.
- The format of the shape files, and `isSurface`:
  [syntaxShapes](https://richefeu.github.io/rockable/syntaxShapes.html).
- The controls of a driven body (velocities, forces, rotations, moments):
  [drivingSystem](https://richefeu.github.io/rockable/drivingSystem.html).
