## 10. Parametric studies with `sweepable`

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

The system of [`helloworldPlus`](../helloworldPlus) is computed eight times,
for four restitution coefficients and two friction coefficients between the
grains and the plane. The eight input files are written by the application
`sweepable` from a short script, and the computations are run in parallel by
`batchable.sh`. The animation shows four of them side by side.

### Notions

- The application **`sweepable`** and its script language (`.foreach`,
  `.read`, `.findLineStartingWith`, `.replaceBy`, `.createFolder`,
  `.saveInFolder`, `.saveCollection`)
  ([sweepable](https://richefeu.github.io/rockable/sweepable.html)).
- The script **`batchable.sh`**, which runs a list of computations in parallel.
- Gathering the results of many computations with a script.

### Files

| File                    | Content                                                         |
|-------------------------|-----------------------------------------------------------------|
| `input.txt`             | The reference conf-file (that of `helloworldPlus`)               |
| `shapes.shp`            | The shapes, shared by all the computations                       |
| `folding_structure.txt` | The script of `sweepable`: the parameters and the folders         |
| `batchable.sh`          | Runs the computations of a list, several at a time               |
| `results.sh`            | Gathers a few results of each computation                        |

### Run it

```sh
sweepable folding_structure.txt              # writes the input files
sh batchable.sh all_inputs.txt rockable 4    # 8 computations, 4 at a time: about 10 seconds
sh results.sh
```

### What happens

**Writing the input files.** `folding_structure.txt` holds two nested loops:

```
.foreach en2Contact 0.01 0.2 0.4 0.6
  .foreach muContact 0.1 0.5
    .read input.txt
    .findLineStartingWith "en2Contact 0 1"
    .replaceBy "en2Contact 0 1 ${en2Contact}"
    ...
    .createFolder "simulations/en2_${en2Contact}_mu_${muContact}"
    .saveInFolder
  .endforeach
.endforeach
.saveCollection "all_inputs.txt"
```

For each pair of values, `input.txt` is read, the lines `en2Contact 0 1` and
`muContact 0 1` (the contacts between the plane, group 0, and the grains,
group 1) are replaced, the line `shapeFile` is changed to point to the shared
`shapes.shp` two folders up, a comment is appended, and the file is saved in
its own folder. `sweepable folding_structure.txt` writes

```
simulations/en2_0.01_mu_0.1/input.txt
simulations/en2_0.01_mu_0.5/input.txt
...
simulations/en2_0.6_mu_0.5/input.txt
all_inputs.txt
```

and `all_inputs.txt` lists the eight input files. Run `sweepable` without
argument to get a fresh copy of `batchable.sh`.

**Running them.** `sh batchable.sh all_inputs.txt rockable 4` goes into the
folder of each file of the list, runs `rockable input.txt` there, with at most
4 computations at the same time, and writes the output of each one in its
`run.log`. The command can be any executable that takes the input file as its
last argument (`rockable -j 2` for instance, between quotes).

**Gathering the results.** `results.sh` reads, for each computation of the
list, the values of `en2` and `mu` in its last conf-file, the final position
`x` of the grain of rice, and the kinetic energy at the end of
`kineticEnergy.txt`:

```
# en2 mu x_rice kinetic_energy
0.01 0.1 -0.509 0.141
0.01 0.5 -0.636 0.37
...
0.6 0.1 -0.539 12.2
0.6 0.5 -0.604 4.66
```

With `mu = 0.5`, the grain of rice ends 6 to 13 cm further to the left than
with `mu = 0.1`. The kinetic energy left at `t = 1 s` is the largest for the
most elastic contacts (`en2 = 0.6`), but it does not vary steadily with `en2`:
each computation follows its own sequence of impacts.

### What to look at

- `see` in any folder of `simulations`, for instance
  `cd simulations/en2_0.6_mu_0.1 && see`.
- `run.log` in each folder, if a computation stopped.
- The comment appended at the end of each input file, which tells with which
  parameters it was written.

### Exercises

1. Add a third value of friction, `.foreach muContact 0.1 0.5 0.9`, and run the
   sweep again: twelve computations.
2. Sweep the time step instead (`.foreach dt 1e-6 2e-6 5e-6`), with
   `.findLineStartingWith dt`: does the result change?
3. In `input.txt`, set `interVerlet 0.01` (the list of neighbours is rebuilt
   ten times less often), then write and run the sweep again. In the case
   `en2 = 0.4, mu = 0.5`, the grain of rice is thrown up at 18 m/s at
   `t = 0.39 s`: its tip, which turns fast, went 1.6 cm into the plane before
   the contact was found, and the repulsive force created energy. A parametric
   study is a good way to catch such numerical artefacts: one case out of
   eight is far from the others (look at its kinetic energy).
4. Plot the final position of the rice against `en2`, for each friction:

   ```sh
   sh results.sh > results.txt
   gnuplot -p -e "plot 'results.txt' using 1:(\$2==0.1?\$3:1/0) title 'mu = 0.1', '' using 1:(\$2==0.5?\$3:1/0) title 'mu = 0.5'"
   ```

### Going further

- The whole script language of `sweepable`, and its C++ interface:
  [sweepable](https://richefeu.github.io/rockable/sweepable.html).
- On a cluster, `batchable.sh` runs as it is; with a job scheduler, the list
  `all_inputs.txt` can also feed a job array.
