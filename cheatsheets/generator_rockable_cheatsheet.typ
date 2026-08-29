#import "common.typ": *

#show: rockable-sheet.with(
  title: [`generator`, the pre-processor that writes `Rockable` input files],
)

= Principle

#concept-block(body: [
  #align(center, command("generator <command-file>"))

  `generator` reads a text script and *writes files*: shape libraries, input
  files, and the particle lists of a packing. It has no notion of a simulation:
  it is a small language for emitting lines, computing numbers, and placing
  bodies.

  A script opens an output file, prints the header of a conf-file, then lets a
  packing command emit the particle lines. `close` comes back and fills in the
  particle count.

  #note[Leave at least one blank line at the end of the command file, otherwise
  the last command may be missed.]
])

= Output and text

#concept-block(body: [
  #kwtable(
    command("open <filename>"), [Open a file for output. Every `print`,
      `compute` and packing command writes there until `close`. With no file
      open, the output goes to the console.],
    command("close"), [Close the current file. If `Particles:placeHolder` was
      used, the final particle count is written back at that position first.],
    command("print <text>"), [Write `text` (leading and trailing spaces
      trimmed) followed by a newline. This is how raw conf-file lines are
      emitted.],
    command("print> <text>"), [Same, *without* the newline: the next output
      continues on the same line.],
  )

  #inline("Arithmetic")
  #kwtable(
    command("compute <label> <expr>"), [Evaluate `expr` and write
      `label result`. `compute Particles 5^3 + 6` gives `Particles 131`.],
    command("<compute <expr>"), [Evaluate `expr` and write only the result,
      preceded by a space. Meant to follow a `print>`.],
  )

  ```txt
  print> knContact 0 0
  <compute 500000 * 2
  ```
  produces `knContact 0 0 1000000`.
])

= Counting the particles

#concept-block(body: [
  - #command("Particles:placeHolder") ~Reserve the spot where the total particle
    count will be written when `close` is called. Every packing command keeps
    the count up to date.
  - #command("incrementNoParticles <n>") ~Increase the tracked count by hand,
    for the particles you emit with `print` rather than with a packing command.
])

= Transformations

#concept-block(body: [
  A global transformation is carried along the script and applied by the packing
  commands when they place the particles.

  - #command("tranformation:add_translation <x> <y> <z>")
  - #command("tranformation:add_rotation <x> <y> <z> <angle_deg>") ~Rotation of
    `angle_deg` degrees around the axis $(x, y, z)$.
  - #command("tranformation:reset") ~Back to identity.
  - #command("tranformation:reset_translation") /
    #command("tranformation:reset_rotation") ~Reset only one part.

  #note[Mind the spelling: the keyword really is `tranformation:`, without the
  first `s`. `transformation:...` is silently ignored (with an "unknown token"
  message).]
])

= Generating shapes

#concept-block(body: [
  These write a shape definition (a `< ... >` block) into the current file. The
  first argument is the shape `name` used later in the particle lines, and
  `radius` is the Minkowski radius.

  - #command("generateShape:sphere <name> <radius>")
  - #command("generateShape:cube <name> <radius> <sideSize>")
  - #command("generateShape:cuboid <name> <radius> <sx> <sy> <sz>")
  - #command("generateShape:pyramid3 <name> <radius> <sideSize>") ~Triangular
    pyramid.
  - #command("generateShape:rhombicuboctahedron <name> <radius> <sx> <sy> <sz>")
  - #command("generateShape:thin_cylinder <name> <Rin> <Rout> <H> <nbSectors>")
    ~Thin cylindrical shell.
  - #command("generateShape:cuboid_container <name> <radius> <sx> <sy> <sz> <hasTop> <hasBottom>")
    ~Open or closed box; the two flags are `0` or `1`.
  - #command("generateShape:rectangle_xz <name> <radius> <side_x> <side_z>")
  - #command("generateShape:rectangle_xy <name> <radius> <side_x> <side_y>")
  - #command("generateShape:xyz_walls <sx> <sy> <sz> <Rw>") ~Three flat wall
    shapes, one per axis.
])

= Generating packings

#concept-block(body: [
  - #command("addParticle <name> <group> <cluster> <homothety> <px> <py> <pz> <qs> <qx> <qy> <qz>")
    ~One particle, at the given position and quaternion orientation.
  - #command("generatePacking:wallBox <group> <LX> <LY> <LZ> <Rw>") ~Six flat
    wall particles forming a closed box, all in `group`.
  - #command("generatePacking:grid <name> <ox> <oy> <oz> <bx> <by> <bz> <nx> <ny> <nz> <group> <cluster> <homothety> <randQ>")
    ~Fill a box of size $(b_x, b_y, b_z)$ from the origin $(o_x, o_y, o_z)$ with
    a regular $n_x times n_y times n_z$ grid. `randQ` set to `1` gives random
    orientations.
  - #command("generatePacking:grid_clust <...>") ~Same arguments, but each
    particle gets its own incrementing cluster index.
  - #command("generatePacking:RandomClosePacking <...>") ~See below.
])

= `RandomClosePacking`

#concept-block(body: [
  #command("generatePacking:RandomClosePacking <name> <boxShape> <direction> <bx> <by> <bz> <xOBB> <yOBB> <zOBB> <hmin> <hmax> <nbTarget> <solidFraction> <group> <cluster>")

  #kwtable(
    command("boxShape"), [Container geometry: `CUBOID` or `CYLINDER`.],
    command("direction"), [Growth direction: `X`, `Y` or `Z`.],
    command("bx by bz"), [Size of the container.],
    command("xOBB yOBB zOBB"), [Bounding-box dimensions of the shape being
      packed.],
    command("hmin hmax"), [Homothety range, i.e. the size scatter.],
    command("nbTarget"), [Target number of particles.],
    command("solidFraction"), [Target solid fraction, in $[0, 1]$.],
    command("group cluster"), [Group and cluster given to every particle.],
  )

  #note[Any `boxShape` other than `CYLINDER` falls back to `CUBOID`. The
  homothety of each placed particle is recovered from its bounding box, so
  `xOBB`, `yOBB` and `zOBB` must describe the shape at homothety 1.]
])

= A complete script

#concept-block(body: [
  ```txt
  # --- shapes ---
  open shapes.txt
  generateShape:cuboid Cuboid 0.0005 0.03 0.01 0.01
  generateShape:xyz_walls 0.25 0.25 0.25 0.0005
  close

  # --- input file ---
  open input.txt
  print Rockable 21-08-2022
  print t 0
  print tmax 2.0
  print dt 5e-6
  print interVerlet 0.001
  print interConf 0.01
  print
  print forceLaw Default
  print Integrator Beeman
  print DVerlet 0.005
  print dVerlet 0.0025
  print density 0 2700
  print density 1 2700
  print knContact 0 1 1e6
  print en2Contact 0 1 0.02
  print ktContact 0 1 1e6
  print muContact 0 1 0.6
  print krContact 0 1 0.0
  print murContact 0 1 0.0
  print iconf 0
  print shapeFile shapes.txt
  print nDriven 0
  Particles:placeHolder

  generatePacking:RandomClosePacking
  Cuboid
  CUBOID Y
  0.2 0.2 0.3
  0.03 0.01 0.01
  1.0 1.0
  2000
  0.9
  0
  0

  close
  ```

  #note[The arguments of a command may be spread over several lines, as above:
  the parser reads tokens, not lines.]
])

= Recipes

#concept-block(body: [
  #inline("A regular grid inside a box of walls")
  ```txt
  open shapes.txt
  generateShape:sphere Ball 0.005
  generateShape:rectangle_xz Wall 0.001 0.2 0.2
  close

  open input.txt
  print Rockable 21-08-2022
  print shapeFile shapes.txt
  print nDriven 6
  Particles:placeHolder
  generatePacking:wallBox 0 0.2 0.2 0.2 0.001
  generatePacking:grid Ball 0.01 0.01 0.01
    0.18 0.18 0.18  6 6 6  1 0 1.0 1
  close
  ```
  `wallBox` emits six particles, hence `nDriven 6`; they must come *first* in the
  list, so the packing command follows.

  #inline("Placing an inclined plane")
  ```txt
  tranformation:reset
  tranformation:add_rotation 0 0 1 -30
  tranformation:add_translation 0 0.1 0
  generatePacking:grid Ball 0 0 0  0.1 0.1 0.1  3 3 3  1 0 1.0 0
  tranformation:reset
  ```
  The transformation stays in effect until it is reset, so a script that mixes
  placed and unplaced groups should reset it explicitly.
])

= Where `generator` fits

#concept-block(body: [
  ```txt
  generator script.txt      -> shapes.txt + input.txt
  shapeSurvey shapes.txt    -> check and pre-compute the shapes
  rockable input.txt -j 8   -> conf0, conf1, ...
  see                       -> look at the result
  ```

  #note[The other packers of #file("prepro/genesis") (`SpherePacker`,
  `TubePacker`) are standalone tools with their own input formats; `generator`
  is the one that writes complete Rockable input files.]
])
