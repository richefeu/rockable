# Tests

This folder holds the tests run by `ctest`: two unit tests of the r-shape skin
mesh, and non-regression tests that run short simulations and compare the final
configuration with a stored reference.

## Running the tests locally

Tests are not built by default (`ROCKABLE_USE_TESTING` is `OFF`). From the root
of the repository:

```sh
cmake -S . -B BUILD -DCMAKE_BUILD_TYPE=Release -DROCKABLE_USE_TESTING=ON
cmake --build BUILD -j
cd BUILD && ctest --output-on-failure
```

Useful variants:

```sh
ctest -N                              # list the tests without running them
ctest -j 6 --output-on-failure        # run them in parallel
ctest -R 518_poly --output-on-failure # run the tests whose name matches a regex
ctest --rerun-failed --output-on-failure
ctest -V -R 50_prismes                # full output, including the deviation found
```

Only the `rockable` executable and the two `test_Shape*` executables are needed,
so the build can be restricted to them when the GUI or prepro dependencies are
missing:

```sh
cmake -S . -B BUILD -DCMAKE_BUILD_TYPE=Release -DROCKABLE_USE_TESTING=ON \
      -DROCKABLE_COMPILE_SEE=OFF -DROCKABLE_COMPILE_PREPRO=OFF \
      -DROCKABLE_COMPILE_POSTPRO=OFF -DROCKABLE_COMPILE_CONF2VTK=OFF
cmake --build BUILD -j --target rockable test_ShapeSDF test_ShapeMesh
```

## Folder layout

```
test/
├── CMakeLists.txt    declares the tests (included when ROCKABLE_USE_TESTING=ON)
├── input/<case>/     input conf-files, shape library, drivingSystem.txt
├── regression/       reference conf-files, one per regression test
└── ShapeSDF/         unit tests of ShapeSDF and ShapeMesh
```

The same references serve every build type (Release, Debug, …) and platform.

## How a regression test works

A regression test is declared in `CMakeLists.txt` by

```cmake
add_regression_test(<name> <case> <conf> <result> <reference>)
```

At configure time, `input/<case>/` is copied into `BUILD/test/<name>/input/`,
and a file `BUILD/test/<name>/redirection.txt` is written:

```
redirection
input/<case>/
<conf>
```

The keyword `redirection` makes `Rockable::loadConf` read the second line as the
path prefix of the input files (shape library, `drivingSystem.txt`) and load the
conf-file named on the third line.

Each regression test is then two ctest tests, both run in `BUILD/test/<name>/`:

1. `<name>_clean` runs `rockable -c`, which deletes the conf-files of a previous
   run. It is a fixture of the next test, so a simulation that stops early can
   never be compared with a stale result.
2. `<name>` runs

   ```sh
   rockable redirection.txt -n <result> -r test/regression/<reference>
   ```

   The simulation runs until `tmax`, then `compareConf` (in `src/Apps/run.cpp`)
   compares the conf-file `<result>` with the reference. If they differ, it
   logs the first rows that differ, and `rockable` exits with a non-zero
   status, which fails the test.

Since each test has its own folder, tests can run in parallel.

### The comparison

`compareConf` reads both files from the `Particles` line on (the header holds
paths and options that may legitimately differ) and cuts them into the sections
`Particles`, `Interactions` and `Interfaces`.

- Section headers, hence the numbers of particles, interactions and interfaces,
  must be identical, and so must non-numeric tokens (shape names).
- Interactions are saved in the order of their addresses in memory, which is not
  reproducible, so they are matched by their key `i j type isub jsub`.
- Two numbers `a` (new) and `b` (reference) in column `j` of a section match when

  ```
  |a - b| <= tolerance * max(|a|, |b|, S_j)
  ```

  where `S_j` is the largest magnitude in column `j` of the reference section.

`S_j` matters for quantities that are the small difference of large ones, such
as the acceleration of a grain at equilibrium: its round-off error scales with
the contact forces, not with the resultant. The tolerance is `1e-5` by default
and can be changed with `-t`. With `ctest -V`, each test prints the largest
relative deviation it found.

Floating-point results are not bit-for-bit reproducible across compilers,
options and processors: on arm64, GCC fuses `a*b + c` into a single FMA
instruction by default, which changes the last digits. The current references
pass with GCC on arm64 with and without `-ffp-contract=off`, and in Debug, with
a largest deviation of 2.6e-6.

## The tests

| Test                                                  | What it runs                                                          |
|-------------------------------------------------------|-----------------------------------------------------------------------|
| `ShapeSDF`                                            | Unit test of the r-shape signed distance function                     |
| `ShapeMesh`                                           | Unit test of the skin mesher and of its coplanar simplification       |
| `2_prismes_gravity`                                   | 2 prisms falling on a plane, 1.5 s (150 000 steps)                    |
| `50_prismes_gravity`                                  | Heap of 50 prisms at rest in a box, 227 contacts, 0.01 s (1000 steps) |
| `518_poly_bruteforce_default_noPeriodicity`           | 518 polyhedra, `Default` law, 11 steps                                |
| `512_particles_bruteforce_LawDefault_noPeriodicity`   | 512 spheres, `Default` law, 11 steps                                  |
| `512_particles_bruteforce_LawAvalanche_noPeriodicity` | 512 spheres, `Avalanche` law, 11 steps                                |

The unit tests print one `PASS`/`FAIL` line per check, with the measured value
and the tolerance, and return a non-zero status if any check fails.

## Writing a regression test

A regression test is only useful if two correct builds agree on its result. Two
kinds of simulations do not:

- **Chaotic ones.** Grains falling and colliding amplify round-off: 50 prisms
  dropped in a box end up in different heaps after about 0.5 s, depending on
  the compiler flags. This is why `50_prismes_gravity` starts from an
  already settled heap and runs only 1000 steps.
- **Laws with a branch decided by round-off.** On a packing at rest, the
  `Avalanche` law switches between loading and unloading on the sign of the
  change of `dn`, which is round-off: a 518-polyhedra test with this law
  diverged by 4% after 11 steps, and was removed.

To check a new test, run it with two builds, one of them configured with
`-DCMAKE_CXX_FLAGS=-ffp-contract=off` (or in Debug), and look at the largest
deviation printed by `ctest -V`. It should stay well below the tolerance.

## Updating a reference

When a change in the code is *meant* to change the results, regenerate the
reference, check the differences, and commit the new file. For instance:

```sh
cd BUILD/test/518_poly_bruteforce_default_noPeriodicity
../../rockable -c
../../rockable redirection.txt
cp conf3001 ../../../test/regression/regression_518_poly_default_noperiod.conf
```

The result file (`conf3001` here) is the fourth argument of
`add_regression_test`.

## Continuous integration

`.github/workflows/run_test.yml` builds Rockable on `ubuntu-latest` with
`ROCKABLE_USE_TESTING=ON` and coverage flags at each push and pull request on
`main`, runs `ctest --no-tests=error` (a configuration without tests fails
instead of passing silently), and uploads a gcovr coverage report.
