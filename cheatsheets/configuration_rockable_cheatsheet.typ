#import "common.typ": *

#show: rockable-sheet.with(
  title: [Configuration file for `Rockable` (conf-files)],
)

= Anatomy of a conf-file

#concept-block(body: [
  A conf-file is a flat list of #command("keyword <values>") entries, parsed in
  the order they appear. Unknown keywords raise a warning and are skipped, so a
  file may be written in any order --- but keywords that *use* data must come
  after the ones that *define* it (#command("shapeFile") before
  #command("Particles"), #command("Particles") before #command("Interactions"),
  #command("Interactions") before #command("Interfaces")).

  #inline("Mandatory first line")
  ```txt
  Rockable 21-08-2022
  ```
  The date is the format version. A different one only triggers a warning.

  #inline("Comments")
  A token starting with `#`, `/` or `!` discards the rest of the line. The
  keyword #command("EOF") stops the parsing before the end of the file.

  #inline("Arithmetic expressions")
  The value of #command("dt") may be given as an expression between `$`:
  ```txt
  dt $2 * 0.5e-6$
  ```

  #inline("Files read next to the conf-file")
  - #file("shapes.txt") ~The shape library (name set by #command("shapeFile")).
  - #file("drivingSystem.txt") ~Boundary driving and servo-controllers.
  - #file("dataExtractors.txt") ~On-the-fly measurement probes.

  #note[The two last files are *not* saved into the conf-files: they are re-read
  from the running folder at each restart.]
])

= Timing

#concept-block(body: [
  - #command("t <value>") ~Current time.
  - #command("tmax <value>") ~Time at which the computation stops.
  - #command("dt <value>") ~Time step increment.
  - #command("interConf <value>") ~Elapsed time between two conf-file dumps.
  - #command("iconf <int>") ~Number of the current configuration. It is used to
    name the next dumps (`conf1`, `conf2`, ...) and is incremented at each dump.
  - #command("precision <int>") ~Number of digits used when writing the
    conf-files (default is the `toofus` CommBox setting).
])

= Neighbor list (NL)

#concept-block(body: [
  - #command("interVerlet <value>") ~Elapsed time between two rebuilds of the
    neighbor list.
  - #command("DVerlet <value>") ~Distance used to decide whether two
    sphero-polyhedra are neighbors. It is added to the Oriented Bounding Boxes
    (OBBs) before testing their overlap; in other words half of this length is
    added on each side of the OBBs.
  - #command("dVerlet <value>") ~Distance used to decide whether two
    sub-elements (sphere for a vertex, tube for an edge, thick polygon for a
    face) of two sphero-polyhedra are neighbors.
  - #command("dynamicUpdateNL <0|1>") ~When set to 1, the list is also rebuilt
    as soon as the largest displacement since the last update exceeds
    #command("dispUpdateNL"), or the largest rotation exceeds
    #command("angleUpdateNL"). These extra updates do not replace the regular
    ones every #command("interVerlet").
    - #command("dispUpdateNL <value>") ~Displacement threshold.
    - #command("angleUpdateNL <value>") ~Rotation threshold, *in degrees*.
  - #command("preventCrossingLength <value>") ~When positive, stiffens the
    normal repulsion at large overlap to prevent bodies from crossing each
    other. Only written in the conf-file when non-zero.
])

= Neighbor-list and contact strategies

#concept-block(body: [
  - #command("UpdateNL <option>") ~Strategy used to rebuild the neighbor list.
    `<option>` is `bruteForce` (default) or `linkCells`. With `linkCells`:
    - #command("cellMinSizes <x> <y> <z>") ~Minimum cell size in each direction.
    - #command("boxForLinkCellsOpt <0|1>") ~Whether the first driven bodies
      belong to the overall bounding box that gets split into cells.
  - #command("AddOrRemoveInteractions <option>") ~Method used to add/remove the
    sub-interactions of a pair of bodies. `<option>` is `bruteForce` (default)
    or `OBBtree`. The best strategy depends on the complexity of the shapes:
    `OBBtree` pays off for shapes with many vertices, edges and faces.
  - #command("ContactPartnership <model>") ~Distributes the interaction
    parameters over the sub-contacts of a same pair of bodies. `<model>` is
    `None` (default), `NumberWeight`, `OverlapWeight` or `SurfaceWeight`.
  - #command("parallel_mode <mode>") ~OpenMP force-accumulation strategy:
    `DefaultParallelMode` (or `Default`), `InteractionBuffer`, `CellMutex`,
    `WaveMethod`, `WaveMethodBlock`. The number of threads is a command-line
    argument (`rockable -j <n>`), not a conf-file entry.
])

= Shapes, groups and particles

#concept-block(body: [
  - #command("shapeFile <path>") ~Path to the shape library. A file of the same
    name found in the running folder takes precedence. The library is not
    re-read when it is already loaded.
  - #command("density <group> <value>") ~Density (kg/m³) of the particles of a
    given group number. Masses and inertias are recomputed from the shape volume
    and the homothety when the particles are read, so #command("density") must
    appear *before* #command("Particles").
  - #command("nDriven <int>") ~The first `nDriven` particles of the list are
    driven bodies (boundaries). Unless #file("drivingSystem.txt") says
    otherwise, they simply do not move.
  - #command("gravity <vec3>") ~Gravity acceleration vector.
  - #command("Particles <nb>") ~Followed by `nb` lines, one per particle:
    #command("<shapeName> <group> <cluster> <homothety> <pos> <vel> <acc> <Q> <vrot> <arot>")
    with `<pos>`, `<vel>`, `<acc>`, `<vrot>`, `<arot>` three reals each and `<Q>`
    a quaternion given as `w x y z`.

  #inline("Editing the particle list by hand")
  Inside the particle list only, a line starting with `#` is a comment, and a
  line starting with `!` is a *disabled* particle: it is skipped _and_ the
  expected count is decremented, so the number after #command("Particles") does
  not need to be edited.
])

= Interactions and interfaces

#concept-block(body: [
  - #command("Interactions <nb>") ~Followed by `nb` lines:
    #command("<i> <j> <type> <isub> <jsub> <n> <dn> <pos> <vel> <fn> <ft> <mom> <damp>")

    `<type>` is `0` vertex-vertex, `1` vertex-edge, `2` vertex-face, `3`
    edge-edge. `<isub>`/`<jsub>` are the sub-element indices in each shape.
    `<n>` is the unit normal oriented from `<j>` to `<i>`, `<dn>` the normal
    distance (negative when overlapping), `<vel>` the velocity of `<j>` relative
    to `<i>` at the contact point, `<fn>` a scalar and `<ft>`, `<mom>` vectors.

  - #command("Interfaces <nb>") ~Glued interfaces between two bodies. One
    interface per line:
    #command("<i> <j> <nbBonds> <dn0>") then, only if
    #command("ParamsInInterfaces") is 1,
    #command("<kn> <kt> <kr> <fn0> <ft0> <mom0> <power> <Gc>"), then `nbBonds`
    triplets #command("<type> <isub> <jsub>") identifying the bonded
    sub-interactions. A bond that cannot be matched with an existing interaction
    raises a warning and the whole interface is dropped.

  - #command("ParamsInInterfaces <0|1>") ~When 1, each interface carries its own
    parameters (written and read on the interface line) instead of taking them
    from the group-pair table. Required by
    #command("setVariableStickParams") and #command("setStiffnessRatioInterfaces").

  - #command("glue_with_walls <yes|no>") ~Also glue the free bodies to the
    driven ones when the sticking prepro-commands are executed.

  - #command("initSpringJoint <i> <ipos0> <j> <jpos0> <stiffness>")
    ~Adds a linear spring joint between bodies `i` and `j`, anchored at
    `ipos0` and `jpos0` given in the *body frames*.
])

= Periodic cell #periodic-flag

#concept-block(body: [
  - #command("usePeriodicCell <0|1>") ~Activate the tri-periodic cell.
  - #command("h <mat9>") ~Cell matrix (9 reals, row by row); its columns are
    the three cell vectors.
  - #command("vh <mat9>") / #command("ah <mat9>") ~Velocity and acceleration of
    the cell matrix.
  - #command("mh <value>") ~Mass ratio used to build the inertia of the cell.
  - #command("dh <value>") ~Numerical damping of the cell degrees of freedom.
  - #command("cellVelocityCorrection <0|1>") ~Remove the mean velocity drift.
  - #command("cellMomentumCorrection <0|1>") ~Remove the total momentum drift.
  - #command("useKineticStress <0|1>") ~Include the kinetic (fluctuation) term
    in the stress used by the cell servo-control.

  #note[The loading of the cell itself is set in #file("drivingSystem.txt") with
  the #command("PeriodicLoading") keyword. Kinematics are stored in *real*
  coordinates in the conf-files and converted to reduced coordinates internally.]
])

= Deformable particles #soft-flag

#concept-block(body: [
  - #command("useSoftParticles <Young> <Poisson>") ~Enable the homogeneous
    straining of the particles and set the compliance from the Young modulus and
    the Poisson ratio.
])

= What a saved conf-file contains

#concept-block(body: [
  Dumps written by `rockable` are themselves valid input files. They always hold
  the timing, the neighbor-list settings, the strategy option names, the group
  properties, every *defined* interaction parameter, then
  #command("Particles"), #command("Interactions") and #command("Interfaces").

  They never hold the pre-processing commands, the data extractors, nor the
  driving system: those live in the initial input file and in the companion
  files, and are re-applied from there.
])

= A minimal input file

#concept-block(body: [
  ```txt
  Rockable 21-08-2022

  t 0
  tmax 1.0
  dt $2 * 0.5e-6$
  interVerlet 0.01
  interConf 0.05

  DVerlet 0.08
  dVerlet 0.02

  density 0 2700
  density 1 2700
  gravity 0 -9.81 0

  forceLaw Default
  Integrator velocityVerlet

  knContact 0 1 1e7
  en2Contact 0 1 0.02
  ktContact 0 1 1e7
  muContact 0 1 0.9
  krContact 0 1 1e7
  murContact 0 1 0.0

  knContact 1 1 1e7
  en2Contact 1 1 0.2
  ktContact 1 1 1e7
  muContact 1 1 0.9
  krContact 1 1 1e7
  murContact 1 1 0.0

  iconf 0
  nDriven 1
  shapeFile shapes.shp
  Particles 3
  #name g c h pos vel acc Q vrot arot
  Plan 0 0 1 0 -.05 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
  Rice 1 0 1 -.5 .5 0 0 0 0 0 0 0 1 0 0 0 0 0 -10 0 0 0
  Ball 1 0 1 .3 .5 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
  ```
  Group `0` holds the driven plane, group `1` the free bodies; every pair of
  groups that can meet must have its parameters defined.
])
