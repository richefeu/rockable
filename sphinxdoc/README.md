# About this documentation

## How it is written

The first version of this documentation was written by hand, by the developers
of Rockable. That worked as long as the code was small enough for one person to
hold it in their head.

It stopped working as the code grew. Force laws, data extractors, servo
controllers and pre-processing commands accumulated, each adding its own
keywords to the conf-file, and the documentation fell behind: by August 2026,
only 41 of the 88 keywords accepted by the parser appeared anywhere in these
pages, two sections consisted of the word `TODO`, and one trailed off on
"Details are listed below" with nothing below.

Since then, the documentation is maintained with the help of **generative AI**,
working from the source code. The point is not to save writing effort. It is
that a machine can be asked to read every `kwMap` entry, every
`computeInteraction`, every `read` method, and to check them one by one against
what these pages claim — which is exactly the kind of exhaustive, tedious
comparison that a human maintainer stops doing once a project passes a certain
size.

Two consequences follow, and both matter.

**The documentation tracks the code.** A page is written from the implementation,
not from the intention. When the two disagree, it is the implementation that is
described.

**The exercise finds bugs.** Reading the source closely enough to document it
turns up things that no one was looking for. The pass of August 2026 produced,
among others:

- `setStiffnessRatioInterfaces` registered its parser entry under the keyword
  `randomlyOrientedVelocities`, silently overwriting it. Neither command did
  what its name said. *Fixed.*
- The documented keyword `UpdateNLStrategy` did not exist; the parser only ever
  registered `UpdateNL`, so the documented spelling was ignored. *Fixed.*
- The normal viscous force was written $\alpha_n \sqrt{2 m_{eff}} v_n$ in three
  places, while the code computes $2 \alpha_n \sqrt{k_n m_{eff}} v_n$.
  *Fixed.*
- `BodyForce ViscousFluid` builds its force from the squared velocity
  components without restoring their sign, so it does not oppose the motion.
  *Documented with a warning; not fixed.*
- `copyParamsToInterfaces inner` does not copy `kr` and `mom0`, while its
  `outer` branch does. *Documented with a warning; not fixed.*
- With `useSoftParticles`, the uniform transformation was not written to the
  conf-files, so such a run could not be restarted exactly from a dump.
  *Fixed* (section `SoftParticles` of the conf-files).

The pass of October 2026, after the port of the soft particles (DCEM), found
the following:

- A computation restarted from a conf-file did not follow the uninterrupted
  one: the computation of the forces done when a conf-file is loaded added one
  more increment to the tangential forces. With `BCM`, the interfaces read from
  a conf-file could no longer break, their area was lost, and so was the
  tangential displacement of their bonds. *Fixed.*
- `TrackDamage` never flushed its file, which lost its end. *Fixed.*
- The spring joints (`initSpringJoint`) were not written in the conf-files, so
  a restarted computation lost them. *Fixed* (section `SpringJoints`).
- The controls of angular velocity (`_zrot_Vel_`...) rotated the body but left
  its angular velocity at zero, so the contacts did not see the motion of a
  rotating drum or blade. *Fixed.*
- `ROCKABLE_ENABLE_BOUNDARY=ON` does not compile (to be adapted from an older
  version of the code).
- `postpro` crashed whatever its command file: the post-processors registered
  themselves with static objects, which the linker dropped from the static
  library. *Fixed* (explicit registration).
- `preventCrossingLength` had no effect: the lines that use it were commented
  out in the force laws. *Fixed* (reactivated, with a capped factor).
- `en2ContactFromViscRate` converted the damping rate into $e_n$ and stored it
  as $e_n^2$, so the damping actually used was smaller than the one given.
  *Fixed.*
- The tangential viscosity of `StickedLinks` was cumulated over time (as in
  `BCM`), and the $d_n^0$ of its inner bonds was taken as zero when the
  parameters come from the group-pair table. *Fixed.*
- The bonds of `StickedLinks` never broke when `mom0*Bond` (or `ft0*Bond`) was
  not defined: the term $0/0$ made the yield function NaN. This was the case of
  the example `stickVerticesInClusters`. *Fixed* (a term with a zero numerator
  contributes nothing).
- The moment threshold of `Avalanche` was undefined when the free body does not
  rotate at all. *Fixed.*
- The numerical damping used the total velocity, not the velocity relative to
  the deformation of a periodic cell. *Fixed.*
- `sawtooth_shaker` drifted when `t_ini` was not a multiple of the period, and
  moved before `t_ini` (which was only a phase); the slope of `ramp` was
  computed once for the whole run. *Fixed.*
- `ParticleStress` built a wrong transverse vector for its columns `a` and `r`
  (*fixed*), and `ClusterGranulo` ignored the values of `SievingSizes`, and
  classified the sub-clusters by their number of particles minus one. *Fixed*:
  the values are the bounds of the classes, in number of particles.
- `conftovtk` copied the file names of the folder into a buffer of 25
  characters. *Fixed.*
- The `generator` commands were spelt `tranformation:...` in the source. *Fixed*
  (the old spelling is still accepted).
- `tess2shape` is not built by CMake, and `src/Makefile` is not maintained.
- Several pages described things the code does not do: a cohesion force in
  `Default`, the integration of the free components of a driven body (they are
  not integrated), gravity on a force-driven wall (it is not applied), the
  coefficients of the `Beeman` and `RungeKutta4` schemes. *Fixed in the pages.*

Findings of this kind are reported, not silently papered over. A page that
described a bug as if it were a feature would be worse than no page at all.

## What this means when you read a page

These pages describe **what the code does**. Where the behaviour is surprising,
incomplete or wrong, there is a `warning` or a `note` saying so, rather than a
tidy description of what it ought to do.

If you find a page that disagrees with the code, the page is the thing to fix,
and the disagreement is worth reporting: it may well be a bug rather than a
typo.

---

# Instructions for updating this documentation with an AI

The rest of this file is addressed to whoever — human or model — takes on the
next documentation pass. Following it keeps the result consistent with what is
already here.

## The one rule

**The source code is the ground truth. The existing documentation is not.**

Never document a keyword, a parameter or an equation from what another page
says about it, from a comment, or from the name of a variable. Open the file
that implements it and read it. Every claim in these pages should be traceable
to a line of C++.

When the code is genuinely ambiguous, say so in the page, or leave the point
out. Do not invent a plausible answer.

## Where the truth lives

| To document | Read |
|---|---|
| conf-file keywords | `Rockable::initParser` in `src/Core/Rockable.cpp` — the `parser.kwMap["..."]` entries |
| what a saved conf-file contains | `Rockable::saveConf` |
| force laws | `src/ForceLaws/ForceLaw_*.cpp`, methods `init()` and `computeInteraction()` |
| interface breakage | `Rockable::check_breakage_of_interfaces` |
| `drivingSystem.txt` | `DrivingSystem::read` in `src/Core/DrivingSystem.cpp` |
| data extractors | `src/DataExtractors/*.cpp`, method `read()` and its `columnDoc` |
| pre-processing commands | `src/PreproCommands/*.cpp`, methods `addCommand()` and `exec()` |
| body forces | `src/BodyForces/*.cpp`, `read()` and `getForceAndMoment()` |
| post-processors | `src/PostProcessors/*.cpp` and `src/Apps/postpro.cpp` |
| shape-file format | `Shape::read` in `src/Core/Shape.cpp` |
| command-line options | `src/Apps/run.cpp`, `see.cpp`, `conftovtk.cpp` |
| build options | `CMakeLists.txt`, the `option(...)` lines |
| registered class names | `Rockable::ExplicitRegistrations` — the factory keys are the C++ class names |

## Finding the gaps

This command lists the keywords the parser accepts that appear nowhere in the
documentation. It is the quickest way to start a pass, and to check it at the
end.

```sh
cd <repo root>
grep -o 'parser\.kwMap\["[^"]*"\]' src/Core/Rockable.cpp \
  | sed 's/.*\["//;s/"\]//' | sort -u > /tmp/kw_code.txt
grep -o 'kwMap\["[^"]*"\]' src/PreproCommands/*.cpp \
  | sed 's/.*\["//;s/"\]//' | sort -u >> /tmp/kw_code.txt
sort -u /tmp/kw_code.txt -o /tmp/kw_code.txt

cat sphinxdoc/source/*.rst | grep -o '``[A-Za-z][A-Za-z0-9_/]*``' \
  | tr -d '`' | sort -u > /tmp/kw_doc.txt

comm -23 /tmp/kw_code.txt /tmp/kw_doc.txt
```

It should print nothing. Other things worth grepping for: `TODO` and `XXX` in
`source/*.rst`, and factory registrations without a matching page.

## Conventions

**Style follows the subject.** Reference material — keywords, options, tool
arguments — is presented as `list-table` or as definition lists, in the manner
of `syntaxConf.rst` and `generator.rst`. Physical models are narrative, with
LaTeX equations, in the manner of `forceLaws.rst` and `integrationSchemes.rst`.
Do not turn one into the other.

**Signatures.** A keyword is introduced with its argument types, in the style
already used throughout:

```rst
``homothetyRange`` (*int*) **ifirst** (*int*) **ilast** (*double*) **hmin**
```

**Unfinished features are documented, with a warning.** A user who enables
`ROCKABLE_ENABLE_BOUNDARY` and sees nothing happen is better served by a page
saying that `Ball::read` and `Cylinder::read` are empty stubs than by silence.
Use `.. warning::` for anything incomplete, untested or surprising, and say
concretely what the code does.

**Cross-reference rather than repeat.** Each page carries a label
(`.. _pageName:`) and is reached with `:ref:`. When a topic already has a page,
point to it instead of restating it — `syntaxConf.rst` points to
`drivingSystem.rst` for exactly this reason.

**Add long pages to a `.. contents::` local toctree**, and register every new
page in the right caption group of `index.rst`, or it will not be reachable.

**Explain the why, not only the what.** "`mh` sets the mass ratio of the cell"
is a restatement of the keyword. "`mh` sets how heavy the cell degrees of
freedom are, and therefore how fast the cell responds; a cell that oscillates
under a constant pressure usually has an `mh` too small" is documentation.

## Verifying

The build must be clean, with **zero warning**, and it must be built **from
scratch**:

```sh
cd sphinxdoc
rm -rf build/html
sphinx-build -b html source build/html
```

The `rm -rf` is not optional. Sphinx only reprocesses the files it thinks have
changed, so an incremental build silently keeps quiet about warnings in the
pages it skipped. Two broken `:ref:` survived several incremental builds during
the August 2026 pass and only appeared on the first clean one.

A `WARNING: label non défini` means a `:ref:` points at a label that does not
exist — usually because the target page has no `.. _label:` at its top.
`WARNING: Le document n'est inclus dans aucune toctree` means a page was not
registered in `index.rst`.

Where a claim can be tested, test it. The example in the top-level `README.md`
was checked by extracting its code blocks and running `rockable` on them, which
is how three errors in it were found. The same applies here: if a page says a
keyword does something, running a two-body case is cheap.

## Committing

Commit in coherent batches — one subject area per commit, not one file — and
write in the message *what was verified*, not only what was added. State
plainly which findings were fixed and which were only documented.

Do not push without being asked.
