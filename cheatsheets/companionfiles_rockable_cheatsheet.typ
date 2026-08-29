#import "common.typ": *

#show: rockable-sheet.with(
  title: [`Rockable` companion files: driving, extractors, prepro and postpro],
  font-size: 5.9pt,
  line-skip: 5pt,
)

= The files around a simulation

#concept-block(body: [
  #kwtable(
    file("input.txt"), [The conf-file you write by hand (any name works).],
    file("shapes.txt"), [Shape library, named by #command("shapeFile").],
    file("drivingSystem.txt"), [Driving and servo-controllers. Read at start-up
      then *re-read at every conf dump*: it can be edited while running.],
    file("dataExtractors.txt"), [On-the-fly measurements. Read once, at
      start-up.],
    file("conf0, conf1, ..."), [Configuration dumps, themselves valid inputs.],
    file("extractedDataDoc.txt"), [Auto-generated: describes each extractor
      file, column by column.],
    file("perf.txt"), [`t`, efficiency (steps/s), and profiling data.],
    file("kineticEnergy.txt"), [`t`, translational energy, rotational energy.],
    file("staticBalance.txt"), [`t`, $F_"max" \/ f_n^"max"$,
      $F_"max" \/ macron(f_n)$ over the free bodies --- the equilibrium check.],
    file("checkplots.txt"), [Auto-generated `gnuplot` script plotting the three
      files above: `gnuplot checkplots.txt`.],
  )

  #note[`rockable -c` deletes `conf*`, `perf.txt`, `kineticEnergy.txt`,
  `staticBalance.txt` and `checkplots.txt` from the current folder, then exits.]
])

= `drivingSystem.txt` --- direct controls

#concept-block(body: [
  Only the first #command("nDriven") bodies can be driven; without this file
  they do not move. Lines starting with `/`, `#` or `!` are comments.

  #align(center, command("Control <type> <bodyNumber> <value>"))

  #kwtable(
    command("_x_Vel_ _y_Vel_ _z_Vel_"), [Imposed velocity component.],
    command("_xrot_Vel_ _yrot_Vel_ _zrot_Vel_"), [Imposed angular velocity
      component.],
    command("_x_For_ _y_For_ _z_For_"), [Imposed force component.],
    command("_xrot_Mom_ _yrot_Mom_ _zrot_Mom_"), [Imposed moment component.],
    command("_xyzrot_Vel_"), [Imposed angular velocity *vector*: the value is
      then three reals instead of one.],
    command("_xyzrot_Mom_"), [Imposed moment *vector*, three reals.],
  )

  ```txt
  Control _y_Vel_ 0 -0.01
  Control _x_For_ 1 -250.0
  Control _xyzrot_Vel_ 2 0 0 1.5
  ```
  A component that is not controlled stays free: the body is integrated
  normally along that degree of freedom.
])

= `drivingSystem.txt` --- servo-controllers

#concept-block(body: [
  #align(center, command("Servo <name> <parameters>"))

  A servo re-computes the values of its controls at every step. Only *one*
  servo is active: the last one read wins.

  #inline("The tritri family (six walls)")
  All start with the six body numbers
  #command("<idXmin> <idXmax> <idYmin> <idYmax> <idZmin> <idZmax>"). The `min`
  walls get a velocity control and the `max` walls a force control, recomputed
  from the current wall spacing so that the *stress* stays at its target.

  - #command("tritriIsostaticCompression <pressure>") ~Same pressure on the
    three `max` walls; the `min` walls are held fixed.
  - #command("tritriBiaxialCompression <pressure> <velocity>") ~Lateral pressure
    along $x$ and $z$; the `Ymax` wall moves down at `velocity`.
  - #command("tritriCustom <xminType> <xminValue> ... <zmaxType> <zmaxValue>")
    ~Six pairs, one per wall. `Type` is `0` for a velocity, `1` for a stress
    (converted to a force using the current wall area).
  - #command("tritriLodeAngle <pressure> <LodeAngle> <sigRate>") ~Stress path at
    constant Lode angle, in $[0°, 60°]$, with $sigma$ increasing at `sigRate`.
    #note[Requires $t = 0$ at the start. Not fully tested.]

  #inline("Shakers")
  - #command("shaker <body> <dir> <A> <freq>") ~Sinusoidal motion of amplitude
    `A` along the (normalised) direction `dir`.
  - #command("triangle_shaker <body> <dir> <A> <freq>") ~Same, triangular wave.
  - #command("sawtooth_shaker <body> <dir> <A> <freq> <t_ini>") ~Triangular wave
    with a settable phase origin `t_ini`.

  #inline("Ramp")
  - #command("ramp <type> <body> <valueBegin> <valueEnd> <tBegin> <tEnd>")
    ~Linear ramp of one control between two times; constant outside.
])

= Periodic loading #periodic-flag

#concept-block(body: [
  In #file("drivingSystem.txt"), driving the periodic cell instead of walls:

  - #command("PeriodicLoading IsotropicCompression <pressure>") ~The three
    diagonal components of the cell are stress-driven, the shear components are
    held at zero velocity.
  - #command("PeriodicLoading TriaxialCompression <X|Y|Z> <pressure> <strainRate>")
    ~The named direction is strain-driven at `strainRate`, the two others are
    kept at `pressure`.
  - #command("PeriodicLoading SimpleShearDeformable <XY|XZ|YX|YZ|ZX|ZY> <pressure> <shearRate>")
    ~Shear at constant `shearRate` on the named component, the normal components
    being stress-driven at `pressure`.
])

= `dataExtractors.txt`

#concept-block(body: [
  One extractor per entry, name first, then its parameters. All of them end with
  #command("<filename> <nrec>"), where `nrec` is a number of *time steps*
  between two records. Comments are allowed (`/`, `#`, `!`).

  - #command("MeanVelocity <file> <nrec>") ~`t`, mean of $norm(v)$ over all the
    bodies.
  - #command("TrackBody <ibody> <file> <nrec>") ~`t`, then `pos`, `vel`, `Q`,
    `vrot`, the resultant force (imposed part removed) and the moment.
  - #command("TrackRockfall <ibody> <vStop> <wStop> <file> <nrec>") ~Columns of
    `TrackBody`, and *stops the run* once the body is slower than
    `vStop`/`wStop` or leaves the bounding box.
  - #command("ClusterAABB <icluster> <file> <nrec>") ~`t` and the six bounds of
    the axis-aligned box of one cluster.
  - #command("dnStat <file> <nrec>") ~Normal-distance statistics: min, max,
    mean, means and counts over the negative and positive ones, and the
    position of the minimum.
  - #command("DuoBalance <i> <j> <file> <nrec>") ~Counts of the four
    sub-interaction types between two bodies, and their weighted sum.
  - #command("TrackDamage <file> <nrec>") ~`t` and the damage, i.e. the broken
    interface area over the initial one.

  #note[Extractors are ignored by the interactive tools (`see`, `conftovtk`). A
  #command("DataExtractor ...") line inside a conf-file still works but is
  deprecated: it would not survive the next dump.]
])

= Pre-processing commands

#concept-block(body: [
  They act on a system already set, so they belong *at the end* of the input
  file. They run once and are never written back into the dumps.

  #inline("Creating glued interfaces")
  - #command("stickVerticesInClusters <epsilon>") ~Glue the bodies of a same
    `cluster`, vertex-to-vertex (sphere bonds) only, closer than `epsilon`.
  - #command("stickVerticesInClustersMoments <epsilon>") ~Same, with the moment
    transmission enabled at the bonds.
  - #command("stickClusters <epsilon>") ~Glue bodies of *different* clusters,
    all bond types, closer than `epsilon`.
  - #command("stickBCM <epsilonDist>") ~Build the interfaces for the `BCM` force
    law (bond areas are computed).

  #note[#command("glue_with_walls yes") extends these commands to the driven
  bodies.]

  #inline("Tuning the interfaces")
  - #command("copyParamsToInterfaces <inner|outer>") ~Copy the group-pair
    parameters into each interface. Needs #command("ParamsInInterfaces 1").
  - #command("setStiffnessRatioInterfaces <ratio>") ~Set $k_t \/ k_n$ on the
    interfaces.
  - #command("setVariableStickParams <param> <inner|outer> <lambda> <m> <0|1>")
    ~Draw one interface parameter from a Weibull law of scale `lambda` and
    modulus `m`. The last flag seeds the generator on the clock.

  #inline("Setting the kinematics")
  - #command("setAllVelocities <vec3>") ~Same velocity for every free body.
  - #command("randomlyOrientedVelocities <magnitude>") ~Random direction, given
    magnitude, per body.
  - #command("randomlyOrientedVelocitiesClusters <magnitude> <opt>") ~Same, but
    one direction per cluster.

  #inline("Populating the sample")
  - #command("homothetyRange <ifirst> <ilast> <hmin> <hmax> <0|1>") ~Random
    homothety in $[h_"min", h_"max"]$ for a range of bodies; the last flag seeds
    on the clock.
  - #command("particlesClonage <ifirst> <ilast> <translation>") ~Duplicate a
    range of bodies, translated by the given vector.
])

= Post-processing (`postpro`)

#concept-block(body: [
  #align(center, command("postpro <commandFile>"))

  The command file selects a range of dumps and one post-processor:

  ```txt
  firstConf 0
  lastConf 100
  stepConf 5
  PostProcessor ParticleStress
  Volume 1.0
  ```

  #kwtable(
    command("ParticleStress"), [Per-particle stress tensor. #command("Volume")
      gives the default total volume, and #command("ConfVolumes <nb>") followed
      by `nb` pairs #command("<iconf> <volume>") overrides it conf by conf.],
    command("ClusterGranulo"), [Size distribution of the (possibly broken)
      clusters. #command("SievingSizes <nb>") followed by `nb` sieve sizes.],
  )

  #note[`postpro` loads the dumps in interactive mode: the extractors are not
  run again.]
])
