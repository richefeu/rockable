## 4. Generating shapes and packings with `generator`

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

3375 small cubes, placed on a grid with random orientations, fall into a box
of 25 cm. Nothing is written by hand: the shapes, the walls, the cubes and the
whole input file are produced by the application `generator`, from short
command files.

### Notions

- The application **`generator`** and its commands `generateShape:...`,
  `generatePacking:...` ([generator](https://richefeu.github.io/rockable/generator.html)).
- Its commands `open`, `print`, `close` and `Particles:placeHolder`, which
  write a complete input file.
- Several shape files, or several pieces of an input file, put together.

### Files

| File              | Content                                                          |
|-------------------|------------------------------------------------------------------|
| `genCube.txt`     | Command of `generator` that writes the shape of a cube            |
| `genBox.txt`      | Command that writes the shapes of the walls of a box              |
| `genPacking.txt`  | Commands that place the 6 walls and the 3375 cubes                |
| `input_head.txt`  | The beginning of the input file, without the particles            |
| `genAll.txt`      | All of the above in a single command file                         |

### Run it

Either step by step:

```sh
generator genCube.txt > shapes.txt
generator genBox.txt >> shapes.txt
generator genPacking.txt > sample.txt
sed 's/^Particles XXXX/Particles 3381/' input_head.txt > input.txt
cat sample.txt >> input.txt
rockable input.txt -j 4          # about 15 minutes on one thread, less with -j 4
see
```

or all at once:

```sh
generator genAll.txt             # writes shapes.txt and input.txt
rockable input.txt -j 4
```

### What happens

**The cube.** `genCube.txt` holds a single command:

```
generateShape:cube Cube 0.0005 0.01
```

a cube named `Cube`, of side 1 cm (the outer size, the rounding included),
with a Minkowski radius of 0.5 mm. `generator` writes the shape on its
standard output, with its volume and inertia (`preCompDone y`).

**The walls.** `genBox.txt`,

```
generateShape:xyz_walls   0.25 0.25 0.25   0.0005
```

writes three wall shapes, `x-wall`, `y-wall` and `z-wall`, normal to `x`, `y`
and `z`, for a box of 25 cm, with a radius (half-thickness) of 0.5 mm. `>>`
appends them to `shapes.txt`.

**The packing.** `genPacking.txt` holds two commands, whose arguments are
given on the following lines:

```
generatePacking:wallBox
1                  group of the walls
0.25 0.25 0.25     inner size of the box
0.0005             half-thickness (the radius of the wall shapes)

generatePacking:grid
Cube               shape
0 0 0              corner of the grid
0.25 0.25 0.25     size of the grid
15 15 15           number of particles along x, y, z
0                  group
0                  cluster
1.0                homothety
1                  1: random orientations, 0: all aligned
```

The six walls come first: with `nDriven 6`, they are the driven bodies, and
since nothing drives them, they stay fixed. The cubes, on a grid of
15 × 15 × 15, are 1.67 cm apart. `generator` tells how many bodies each command
added (on the error output, so that it does not end up in the file):

```
@generatePacking:wallBox, Number of added bodies: 6
@generatePacking:grid, Number of added bodies: 3375
```

**The input file.** `input_head.txt` holds the parameters of the computation
and ends with `Particles XXXX`: `sed` puts the number of particles, 3381, in
its place, and the lines of `sample.txt` are appended.

**All at once.** `genAll.txt` does the same in a single file. `open shapes.txt`
sends what follows into `shapes.txt` until `close`; `print` writes a line as it
is; `Particles:placeHolder` writes `Particles` and leaves room for the number,
which `generator` fills in at `close` with the count of the particles added in
between.

**The computation.** There is no `gravity` line: the default gravity,
`(0, -9.81, 0)`, applies. The cubes fall, pile up at the bottom of the box,
and are at rest from about `t = 0.5 s`. With 3375 particles, the computation
takes about 15 minutes on one thread; `-j 4` runs it on 4 threads.

### What to look at

- In `see`, the key `d` hides the walls.
- The kinetic energy, `kineticEnergy.txt`, which drops once the cubes are
  piled up:

  ```sh
  gnuplot -p -e "set logscale y; plot 'kineticEnergy.txt' using 1:2 with lines title 'translation'"
  ```
- `shapes.txt`: the shapes written by `generator`.

### Exercises

1. Replace the cube by an elongated cuboid: `generateShape:cuboid Cube 0.0005
   0.015 0.005 0.005` (the name is kept, so that nothing else has to change).
2. Put fewer cubes, aligned (`5 5 5` and a last `0` in `generatePacking:grid`).
   The number of particles is then 131: change it in the command `sed`, or use
   `genAll.txt`, which counts them.
3. Make the box twice as high: `0.25 0.5 0.25` in `genBox.txt`, and for the
   size of the box and of the grid in `genPacking.txt`.
4. Look at the other commands of `generator` (spheres, polyhedra, random
   close packings, transformations) on the page
   [generator](https://richefeu.github.io/rockable/generator.html).

### Going further

- [`DominoFun`](../DominoFun): a whole system written by `generator`, with
  clusters of glued particles.
- The C++ functions behind `generator`, in `prepro/genesis/cpptools`, can be
  used directly in a program.
