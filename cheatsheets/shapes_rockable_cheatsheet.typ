#import "common.typ": *

#show: rockable-sheet.with(
  title: [Shape files, `shapeSurvey` and the converters],
)

= What a shape is

#concept-block(body: [
  A `Rockable` particle is a *sphero-polyhedron* (R-shape): the Minkowski sum of
  a skeleton --- vertices, edges, faces --- with a ball of radius $R$. There is
  a single radius per shape, so all its edges and corners are rounded the same
  way. The skeleton may be non-convex, may have holes, and may even be an open
  surface.

  A shape-file is a *library*: several shapes in one file, each written between
  a `<` and a `>`. All the shapes a simulation refers to must be in the single
  file given by #command("shapeFile").

  #note[A shape-file of the same name found in the running folder takes
  precedence over the path stored in the conf-file --- handy to override a
  library locally.]
])

= Keywords of a shape block

#concept-block(body: [
  #inline("Identity")
  - #command("name <string>") ~Name, without spaces. This is what the particle
    lines of the conf-file refer to.
  - #command("radius <double>") ~The Minkowski radius $R$.
  - #command("preCompDone <y|n>") ~When `y`, the mass properties are trusted as
    written and not recomputed at start-up.
  - #command("isSurface") ~Flag (no argument): the shape is an open surface, not
    a solid.

  #inline("Skeleton")
  - #command("nv <int>") ~Then that many lines #command("<x> <y> <z>").
  - #command("ne <int>") ~Then that many lines #command("<from> <to>"), vertex
    indices.
  - #command("nf <int>") ~Then that many lines
    #command("<nbVertices> <v1> <v2> ...").

  #inline("Mass properties")
  - #command("volume <double>") ~Volume of the shape.
  - #command("I/m <I1> <I2> <I3>") ~Eigenvalues of the inertia tensor divided by
    the mass. It assumes the vertices are given in the eigen-frame.
  - #command("MCnstep <int>") ~Number of Monte-Carlo samples used when the mass
    properties have to be computed.

  #inline("Bounding box")
  - #command("obb.extent <e1> <e2> <e3>") ~Half-extents along the three axes.
  - #command("obb.e1 <vec3>"), #command("obb.e2 <vec3>"),
    #command("obb.e3 <vec3>") ~The three unit axes.
  - #command("obb.center <vec3>") ~Centre of the OBB relative to the mass
    centre, in the frame of the shape.
  - #command("fibObbOption <int>") ~Strategy used to fit the OBB: `0`
    covariance, `1` minimum volume, `2` axis-aligned, `3` imposed axis.

  #inline("Optional reminders")
  - #command("position <vec3>") / #command("orientation <quat>") ~A position and
    an orientation kept for processing purposes. Usually written by an automated
    pre-computation, rarely set by hand.
  - #command("OBBtreeLevel <int>") ~(deprecated) Number of levels of the
    OBB-tree.
])

= A unit cube, radius 0.1

#concept-block(body: [
  ```txt
  <
  name Cube_r0.1
  radius 0.1
  preCompDone y

  nv 8
  0.4 0.4 -0.4
  -0.4 0.4 -0.4
  -0.4 -0.4 -0.4
  0.4 -0.4 -0.4
  0.4 0.4 0.4
  -0.4 0.4 0.4
  -0.4 -0.4 0.4
  0.4 -0.4 0.4

  ne 12
  0 1
  1 2
  2 3
  3 0
  4 5
  5 6
  6 7
  7 4
  0 4
  1 5
  2 6
  3 7

  nf 6
  4 0 1 2 3
  4 4 5 6 7
  4 0 1 5 4
  4 2 3 7 6
  4 1 2 6 5
  4 0 4 7 3

  obb.extent 0.5 0.5 0.5
  obb.e1 1 0 0
  obb.e2 0 1 0
  obb.e3 0 0 1
  obb.center 0 0 0

  volume 0.975587
  I/m 0.166667 0.166667 0.166667
  >
  ```

  The skeleton is a cube of side $0.8$; dilated by $R = 0.1$ it becomes a
  rounded cube of overall side $1$. That is why the vertices sit at
  $plus.minus 0.4$, and why the volume is slightly below $1$.
])

= `shapeSurvey`

#concept-block(body: [
  #align(center, command("shapeSurvey <shapeFile>"))

  An OpenGL browser for a shape library: check that a shape is correct, compute
  its mass properties, fit its OBB, then save the library back.

  #inline("Browsing")
  - #command("+ / -") ~Next / previous shape. A shape whose `preCompDone` is `n`
    gets its OBB fitted on the fly.
  - #command("h") ~Show the help. #command("q") ~Quit.
  - #command("e") ~Print the OBB extents in the terminal.
  - #command("a / A") ~Decrease / increase the transparency.
  - #command("b") ~Background on/off.
  - #command("w / W") ~Roll the camera around the view axis.

  #inline("Pre-computing")
  - #command("c") ~Compute the mass properties of the current shape and set its
    `preCompDone` to `y`.
  - #command("C") ~Same, for *every* shape of the library that still has `n`.
  - #command("*") ~Reset `preCompDone` of the current shape to `n`, to force a
    recomputation.
  - #command("N / n") ~Multiply / divide by 10 the number of Monte-Carlo steps
    (`MCnstep`), between $10^3$ and $10^8$. More steps means a more accurate
    volume and inertia.
  - #command("o") ~Cycle the OBB fitting strategy: covariance, minimum volume,
    axis-aligned, imposed axis.
  - #command("t") ~Build the OBB-tree; #command("k / K") ~show one level less /
    more of it.
  - #command("d") ~Clean every shape (remove the duplicated entities).

  #inline("Saving")
  - #command("s") ~Save the library back to *the file it was read from*.
  - #command("p") ~Export a sample: one particle per shape, laid out so that the
    library can be opened directly in `see`.

  #note[The usual workflow after a conversion: open the library, press
  #command("C"), then #command("s"). The simulation then starts without
  recomputing anything.]
])

= `stl2shape`

#concept-block(body: [
  #align(center, command("stl2shape -i <file.stl> -r <radius> [options]"))

  Converts a *binary* STL mesh into one shape.

  #kwtable(
    command("-i, --input <string>"), [Input STL file. Required.],
    command("-r, --radius <double>"), [Minkowski radius of the result. Required.],
    command("-s, --scaleFactor <double>"), [Rescale the object.],
    command("-m, --maxLength <double>"), [Maximum length (sieving size) of the
      result; it is rescaled to match.],
    command("-z, --scaleRadius"), [Also rescale the radius along with the
      object.],
    command("-c, --clean"), [Remove the duplicated edges. Almost always worth
      it: an STL repeats every edge.],
  )

  #note[A mesh that is too fine gives a shape with thousands of sub-elements,
  and the contact detection cost follows. Decimate before converting, and
  consider #command("AddOrRemoveInteractions OBBtree").]
])

= `tess2shape`

#concept-block(body: [
  #align(center, command("tess2shape <command-file>"))

  Turns a `neper` tessellation (`.tess`) into a shape library *and* the
  matching particle list --- the usual route to a polycrystal-like sample of
  space-filling grains.

  ```txt
  tessFileName        tessel.tess
  inputFileName       input.txt
  shapesFileName      shapes.txt
  MinkowskiRadius     0.01
  ParticlesGroup      0
  ParticlesCluster    0
  ParticlesHomothety  1.0
  MCnstep             50000
  ```

  #note[The keyword is `MinkowskiRadius`. Some example files in the repository
  spell it `MinskowskiRadius`, which is *not* recognised and silently leaves the
  radius at its default.]

  Generating the tessellation itself:
  ```sh
  neper -T -n 100 -domain "cube(1,1,1)" -o tessel
  ```
])

= `shape2mesh`

#concept-block(body: [
  #align(center, command("shape2mesh <input.shp> [options]"))

  Meshes the *exact* surface of the dilated solid --- rounded edges, corners,
  concavities and open surfaces included --- with analytic vertex normals, so
  even a coarse mesh shades well. The mesh is expressed in the body frame.

  #kwtable(
    command("-f, --format obj|ply"), [Output format (default `obj`).],
    command("-e, --epsilon <value>"), [Absolute sag tolerance, in the length
      unit of the `.shp`. Smaller is finer; the grid step follows
      $h tilde.op sqrt(8 epsilon R)$. Default: 1% of the radius.],
    command("-s, --simplify [deg]"), [Merge the coplanar triangles (default
      angle $0.25°$) and re-triangulate the flat regions from their boundary.
      The mesh stays watertight; curved regions are untouched.],
    command("-r, --rmsh"), [Write a single #file("<stem>.rmsh") holding every
      shape of the library.],
    command("-o, --outdir <dir>"), [Output directory.],
  )

  #inline("The .rmsh companion")
  ```sh
  shape2mesh shapes.txt --rmsh       # -> shapes.rmsh
  shape2mesh shapes.txt --rmsh -s    # smaller, flat faces merged
  ```
  With a #file("shapes.rmsh") next to #file("shapes.txt"), `see` draws each
  particle with its skin mesh instead of the overlapping primitives. Without an
  explicit `-e`, the sag is chosen per shape (0.4% of its bounding-box
  diagonal), so a large flat wall and a small detailed grain both get a sensible
  mesh.

  #note[Simplification typically removes ~60% of the triangles of a rounded
  polyhedron and leaves a sphere untouched.]
])
