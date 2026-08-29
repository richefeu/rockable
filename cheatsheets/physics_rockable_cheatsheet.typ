#import "common.typ": *

#show: rockable-sheet.with(
  title: [Force laws, dissipation and time integration in `Rockable`],
)

= Setting an interaction parameter

#concept-block(body: [
  Every interaction parameter is stored in a symmetric table indexed by the
  *group numbers* of the two bodies:

  #align(center, command("<parameter> <group1> <group2> <value>"))

  ```txt
  knContact 0 2 1e8
  ```
  sets $k_n = 10^8$ between the bodies of group `0` and those of group `2`.
  Setting `0 2` also sets `2 0`. A pair that is never defined keeps a null
  value, so *every pair of groups that can meet must be given its parameters*.

  Only the parameters registered by the selected #command("forceLaw") are
  meaningful; the others are simply never read. Parameters are written back to
  the conf-files, but only for the group pairs that were actually defined.
])

= Choosing the force law

#concept-block(body: [
  #command("forceLaw <name>") ~with `<name>` one of:

  #kwtable(
    command("Default"), [Linear elasticity, normal viscosity, Coulomb friction
      and rolling resistance. The one to use in most cases.],
    command("Avalanche"), [Historical law for rock avalanches (inherited from
      `DEMbox`). Normal force with distinct loading/unloading paths.],
    command("StickedLinks"), [Breakable elastic solid bonds, with a stress-like
      rupture criterion. Introduced for the cemented granular materials.],
    command("BCM"), [Bonded Cell Method: bonds carry a stiffness shared over
      the sub-bonds of the interface, with an energy-based (`Gc`) criterion.],
    command("GeoVisc"), [Variant of `Default` for slow landslides: the friction
      force is *viscous* while still yielding on the Coulomb cone.],
  )

  #note[An unknown name falls back to `Default` with a warning.]
])

= Contact parameters

#concept-block(body: [
  Used by all the force laws (`viscTContact` only by `GeoVisc`).

  #kwtable(
    command("knContact"), [Normal stiffness $k_n$. The elastic normal force is
      $f_n^"elas" = -k_n d_n$ for $d_n <= 0$, zero otherwise.],
    command("ktContact"), [Tangential stiffness $k_t$, used incrementally:
      $Delta f_t = [k_t v_t Delta t]_(plus.minus mu f_n^"elas")$.],
    command("muContact"), [Coulomb friction coefficient $mu$.],
    command("krContact"), [Rolling stiffness $k_r$. The resistant moment is only
      computed when $k_r > 0$.],
    command("murContact"), [Rolling resistance coefficient $mu_r$, used as a
      *length*: the moment is capped at $mu_r ell f_n^"elas"$.],
    command("en2Contact"), [Squared normal restitution coefficient $e_n^2$. The
      viscous force is $f_n^"visc" = 2 alpha_n sqrt(k_n m_"eff") v_n$ with
      $m_"eff" = (m_i m_j)/(m_i + m_j)$ and
      $alpha_n = -ln e_n \/ sqrt(ln^2 e_n + pi^2)$. When one of the bodies is
      driven, $m_"eff"$ is the mass of the free one.],
    command("en2ContactFromViscRate"), [Same target, given the other way round:
      the value is the damping rate $alpha_n$, and
      $e_n^2 = exp(-alpha_n pi \/ sqrt(1 - alpha_n^2))$ is stored.],
    command("viscTContact"), [Tangential viscosity of `GeoVisc`:
      $f_t = "viscT" dot v_t$, still capped by $mu f_n^"elas"$.],
  )

  #note[Viscous damping can make $f_n$ negative; it is clamped to zero by
  default. Objectivity correction of $f_t$ under large rotations is a
  compile-time option #ftcorr-flag.]
])

= Bond parameters

#concept-block(body: [
  A bond is *inner* when the two bodies share the same `cluster` number, and
  *outer* otherwise. Replace `*` below by `Inner` or `Outer`.

  #kwtable(
    command("kn*Bond"), [Normal stiffness of the bond,
      $f_n^"elas" = -k_n (d_n - d_n^0)$ with $d_n^0$ the distance at gluing.],
    command("kt*Bond"), [Tangential stiffness (incremental, no yield).],
    command("kr*Bond"), [Rolling stiffness of the bond.],
    command("en2*Bond"), [Squared normal restitution of the bond.],
    command("fn0*Bond"), [Normal strength $f_n^0$ (tension).],
    command("ft0*Bond"), [Tangential strength $f_t^0$.],
    command("mom0*Bond"), [Moment strength $M_0$.],
    command("pow*Bond"), [Exponent $p$ of the rupture criterion.],
    command("gc*Bond"), [Fracture energy $G_c$ (energy-based criterion, `BCM`).],
  )

  #inline("Rupture criterion of StickedLinks")
  $ phi.alt = (abs(f_t) / f_t^0)^p + (norm(M) / M_0)^p - f_n / f_n^0 - 1 $
  The bond breaks irreversibly as soon as $phi.alt > 0$; the breakage of the
  whole interface is then scheduled.

  #note[`Default`, `Avalanche` and `GeoVisc` register only the contact
  parameters. `StickedLinks` registers the full inner/outer sets;
  `BCM` registers `kn`, `kt`, `kr` and `gc` for the inner bonds.]
])

= Time-dependent parameters (`Tempo`)

#concept-block(body: [
  A #command("Tempo") entry plugs a value onto a time-driven function, evaluated
  at every step. Two profiles are available:

  #kwtable(
    command("Range <t1> <t2> <v1> <v2>"), [`v1` while $t in [t_1, t_2]$,
      `v2` outside.],
    command("Ramp <t1> <t2> <v1> <v2>"), [`v1` before $t_1$, linear ramp from
      `v1` to `v2` between $t_1$ and $t_2$, `v2` after $t_2$.],
  )

  #inline("Three targets")
  - #command("Tempo NDCoeff <cmd> <t1> <t2> <v1> <v2>")
    ~Drives #command("numericalDampingCoeff").
  - #command("Tempo Inter <param> <g1> <g2> <cmd> <t1> <t2> <v1> <v2>")
    ~Drives one interaction parameter of a group pair (both `g1 g2` and `g2 g1`).
  - #command("Tempo Body <prop> <group> <cmd> <t1> <t2> <v1> <v2>")
    ~Drives one body property (e.g. `density`) of a group.

  ```txt
  Tempo Inter muContact 0 1 Ramp 0.0 0.5 0.9 0.2
  Tempo NDCoeff Range 0.0 0.1 0.05 0.0
  ```
])

= Global dissipation

#concept-block(body: [
  #inline("Cundall numerical damping")
  #command("numericalDampingCoeff <value>") ~A purely numerical damping that
  stabilises quasi-static simulations. Both translations and rotations are
  affected:
  ```cpp
  if (force * velocity > 0.0) force *= (1.0 - coeff);
  else                        force *= (1.0 + coeff);
  ```
  Typical values range from `0.001` to `0.1` and above; it has to be calibrated,
  and a large value adds artificial dissipation.

  #inline("Velocity barriers")
  A barrier function caps the velocities without adding a dissipative force
  (an idea of Farhang Radjaï). For each component:
  ```cpp
  ratio = pow(fabs(velocity / Barrier), Exponent);
  force *= (1.0 - ratio) / (1.0 + ratio);
  ```
  - #command("VelocityBarrier <value>") + #command("VelocityBarrierExponent <value>")
  - #command("AngularVelocityBarrier <value>") + #command("AngularVelocityBarrierExponent <value>")

  #note[Useful to deposit particles without polluting the contact law.]
])

= Body forces

#concept-block(body: [
  #command("BodyForce <name> <parameters>") ~At most one body force at a time.

  - #command("ViscousFluid <fluidDensity>") ~Drag of a surrounding viscous
    fluid. The damping is $1/2 C_X rho V^(2\/3)$ with $C_X = 0.47$ (sphere).
  - #command("AttractingPoint <point> <acceleration>") ~Constant acceleration
    towards a fixed point, given as three reals.
  - #command("PreferredDirection <axisBody> <axis> <momentMax>") ~Restoring
    moment bringing the body axis `axisBody` (body frame) onto the fixed
    direction `axis`, with $k_r = 2 M_"max" \/ pi$.

  #note[Gravity is *not* a body force: it is the #command("gravity") keyword.]
])

= Time-integration schemes

#concept-block(body: [
  #command("Integrator <name>")

  #kwtable(
    command("velocityVerlet"), [Default. Second-order, time-reversible,
      energy-conserving. New features land here first --- check the code before
      using another scheme.],
    command("Euler"), [First-order, conditionally stable. Errors accumulate;
      mostly useful as a reference.],
    command("Beeman"), [Second-order, uses the acceleration of the previous
      step (an extra initialisation is done at start-up).],
    command("RungeKutta4"), [Runge-Kutta-Nyström, fourth order. Four force
      evaluations per step, better on stiff systems.],
  )

  For all the schemes, the bodies driven in *force* or *moment* are integrated
  like the free ones, whereas the bodies driven in *velocity* are simply moved
  by $x_(t + Delta t) = x_t + v_t^"imposed" Delta t$ (and likewise for the
  orientation quaternion).
])

= Choosing the time step

#concept-block(body: [
  The critical time step of a single contact is
  $ Delta t_c = pi sqrt(m_"eff" \/ k_n) $
  where $m_"eff"$ is the effective mass of the pair (the mass of the free body
  when the other one is driven) and $k_n$ is the relevant normal stiffness
  (`knContact`, or `kn*Bond` for a bonded pair).

  At start-up, `rockable` prints three ratios $Delta t_c \/ Delta t$ (visible at
  verbose level `info`):
  - *estimated* --- from a single contact between two particles;
  - over *ALL* the interactions of the neighbor list;
  - over the *ACTIVE* interactions only.

  Keep $Delta t$ well below $Delta t_c$; the usual practice is a ratio of a few
  tens. In `see`, the key #command("i") displays the same ratio for the
  configuration being viewed.

  #note[The check is skipped for a pair of groups whose `knContact`,
  `knInnerBond` and `knOuterBond` are all undefined --- which is also the sign
  that a parameter is missing.]
])

= A consistent parameter set

#concept-block(body: [
  A minimal, self-consistent contact block for two groups `0` (walls) and `1`
  (grains):

  ```txt
  forceLaw Default
  Integrator velocityVerlet
  numericalDampingCoeff 0.0

  knContact 0 1 1e7
  ktContact 0 1 1e7
  en2Contact 0 1 0.02
  muContact 0 1 0.9
  krContact 0 1 0.0
  murContact 0 1 0.0

  knContact 1 1 1e7
  ktContact 1 1 1e7
  en2Contact 1 1 0.2
  muContact 1 1 0.9
  krContact 1 1 0.0
  murContact 1 1 0.0
  ```

  Points worth checking:
  - the pair `0 0` may be left undefined only if two walls never meet;
  - $k_t \/ k_n$ is usually taken in $[0.2, 1]$;
  - $k_r = 0$ switches the rolling resistance off entirely;
  - `en2Contact` is $e_n^2$, *not* $e_n$;
  - a bonded simulation also needs its `*InnerBond` / `*OuterBond` set, plus a
    prepro-command to actually create the interfaces (see the companion-files
    sheet).
])
