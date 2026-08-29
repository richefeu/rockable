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
- With `useSoftParticles`, the uniform transformation is not written to the
  conf-files, so such a run cannot be restarted exactly from a dump.
  *Documented with a warning; not fixed.*

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
