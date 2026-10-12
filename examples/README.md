# Rockable examples and tutorials

This folder holds small simulations that show how to use Rockable. Most of them
are meant to be followed in order, as the tutorials of a training day: each one
introduces a few notions, and each folder has a `readme.md` that explains the
files, how to run the example, what to look at, and a few exercises.

The reference documentation (all the keywords, the force laws, the file
formats) is the Sphinx documentation, in `sphinxdoc`, also published at
<https://richefeu.github.io/rockable/>. The one-page sheets of the folder
`cheatsheets` are handy to keep at hand during a tutorial.


## Before starting

Build and install Rockable from the root of the repository, then make the
executables reachable from the shell:

```sh
sh install_rockable.sh
source add_install_to_path.sh
```

The examples then use the commands `rockable` (the computation), `see` (the
viewer), and, for some of them, `generator`, `shapeSurvey`, `stl2shape`,
`postpro` or `sweepable`.

An example can be run directly in its folder:

```sh
cd examples/helloworld
rockable input.txt
see conf10
rockable -c          # deletes the results of the computation
```

The files written by a computation (`conf0`, `conf1`, ..., `kineticEnergy.txt`,
`perf.txt`, ...) are ignored by git (see `examples/.gitignore`), so they cannot
be committed by mistake. `rockable -c` deletes them, except the files written
by the data extractors and `extractedDataDoc.txt`, which have to be deleted by
hand. To keep the folder of the
repository clean during a training day, the folder can also be copied first:
`cp -r examples/helloworld ~/tuto/`.


## Learning path

| #  | Folder                                                 | What it teaches                                                                 | Time    |
|----|--------------------------------------------------------|---------------------------------------------------------------------------------|---------|
| 1  | [`helloworld`](helloworld)                             | The input file, the shape file, running a computation, looking at it with `see` | < 1 s   |
| 1b | [`helloworldPlus`](helloworldPlus)                     | A shape from an STL mesh, the homothety, a spring joint, an expression          | 3 s     |
| 2  | [`driving`](driving)                                   | Driving a body with `drivingSystem.txt`                                          | 1 s     |
| 3  | [`isoComprSpheres`](isoComprSpheres)                   | Preparing a sample: a box of walls, a servo-controlled isotropic compression     | 20 s    |
| 4  | [`prepro_cpptools`](prepro_cpptools)                   | Generating shapes and packings with `generator`                                  | short   |
| 5  | [`DominoFun`](DominoFun)                               | A complete system built with `generator`                                         | TBD     |
| 6  | [`stickVerticesInClusters`](stickVerticesInClusters)   | Breakable bonds: clusters, the `StickedLinks` law                                | 4 s     |
| 7  | [`isoComprSpheresPeriodic`](isoComprSpheresPeriodic)   | Periodic boundary conditions                                                     | 20 s    |
| 8  | [`PerioPoly`](PerioPoly)                               | A periodic sample of polyhedra under shear                                       | TBD     |
| 9  | [`softParticles`](softParticles)                       | Deformable particles (DCEM), bonded cells (`BCM`)                                | 1-5 s   |
| 10 | [`sweepable`](sweepable)                               | Parametric studies: many computations from one script                            | short   |
| 11 | [`complexShape`](complexShape)                         | From an STL mesh to a shape: `stl2shape`, `shapeSurvey`                          | short   |

The times are those of the whole computation on a laptop, with one thread.


## Status of the other folders

These folders are kept for now, but they are **not** tutorials. They will be
reviewed, and possibly removed.

| Folder                        | Status                                                                                                                                                |
|-------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`BigForTests`](BigForTests)  | A large bonded sample (23 114 particles, 11 MB) used to test the performance. Its file `input.chai` is **obsolete**: ChaiScript is no longer supported. |
| [`test_GeoVisc`](test_GeoVisc)| A test of the force law `GeoVisc` (2 758 particles), without explanations. To be reviewed.                                                            |
| [`test_mukesh`](test_mukesh)  | Development tests of the `BCM` law (22 MB, with result files). **Superseded** by `softParticles`; candidate for removal.                               |

In `helloworld`, the file `hello.chai` is obsolete for the same reason as
`input.chai` (ChaiScript is no longer supported).


## How each example is documented

The `readme.md` of each tutorial follows the same plan:

1. **Goal** - what the example shows, in two or three lines.
2. **Notions** - the keywords and tools it introduces, with links to the
   documentation.
3. **Files** - what each file of the folder is for.
4. **Run it** - the commands, and how long it takes.
5. **What happens** - the physics of the computation.
6. **What to look at** - in `see`, in the output files, with the expected
   values when there are some.
7. **Exercises** - small changes to try, for a training day.
8. **Going further** - pointers to the documentation and to other examples.

A new example should follow the same plan, and be added to the tables above.
