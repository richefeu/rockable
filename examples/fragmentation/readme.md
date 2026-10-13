## 13. Fragmentation: a block of bonded cells breaks on impact

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

A cube of 6 × 6 × 6 cubic cells, glued face to face, falls on one of its
corners onto a plate at 10 m/s. The interfaces between the cells break when
the elastic energy they store exceeds the energy needed to open a crack: a cone
is crushed under the corner, fragments are ejected, and the top of the block
stays whole. The damage, the fragments and the energy balance are measured.

### Notions

- The force law **`BCM`** (bonded cell method), whose interfaces break with a
  **fracture energy** `Gc`
  ([forceLaws](https://richefeu.github.io/rockable/forceLaws.html)).
- The pre-processing command **`stickBCM`**, which glues the cells of a same
  cluster ([preProcessing](https://richefeu.github.io/rockable/preProcessing.html)).
- **`generatePacking:grid_clust`** and the **transformations** of `generator`
  ([generator](https://richefeu.github.io/rockable/generator.html)).
- The data extractor **`TrackDamage`**, and the post-processor
  **`ClusterGranulo`**, which counts the fragments
  ([dataExtractors](https://richefeu.github.io/rockable/dataExtractors.html),
  [postProcessing](https://richefeu.github.io/rockable/postProcessing.html)).
- **Soft cells** (exercise 4), and the stability of their update
  ([softParticles](https://richefeu.github.io/rockable/softParticles.html)).

### Files

| File                 | Content                                                        |
|----------------------|----------------------------------------------------------------|
| `input.txt`          | The conf-file: the plate and the 216 cells                      |
| `genShapes.txt`      | Commands of `generator` that write `shapes.txt`                |
| `shapes.txt`         | The shapes `Cell` (a rounded cube of 1 cm) and `Platen`         |
| `genBlock.txt`       | Commands of `generator` that write the lines of the cells       |
| `dataExtractors.txt` | `TrackDamage`: the broken fraction of the interfaces (`damage.txt`) |
| `fragments.txt`      | The command file of `postpro`: the fragments in `conf20`        |
| `energies.sh`        | The energy balance of the impact                               |

### Run it

```sh
rockable input.txt       # about 7 seconds, conf0 to conf20
see
postpro fragments.txt    # writes ClusterGranulo-conf20.txt
sh energies.sh
```

### What happens

**The block.** `genShapes.txt` makes a cube of 1 cm with a Minkowski radius of
1 mm (`Cell`), and a plate (`Platen`). `genBlock.txt` places 6 × 6 × 6 cells on
a grid, turned so that a corner points down:

```
transformation:add_translation 0 0.056 0
transformation:add_rotation -1 0 1 54.7356
generatePacking:grid_clust Cell
...
```

The rotation of 54.74° about the axis `(-1, 0, 1)` brings the diagonal of the
cube, `(1, 1, 1)`, along `y`. `grid_clust` gives all the cells the **same
cluster** (0), while `grid` gives each cell its own: `stickBCM` only glues the
cells of a same cluster.

**The glue.** At the end of `input.txt`, `setAllVelocities 0 -10 0` launches the
cells, and `stickBCM 1e-4` glues them. Two neighbouring cells face each other
with four pairs of vertices whose spheres touch (they are less than 1e-4 m
apart): these four **bonds** make one **interface**, whose area `A` is
computed from the positions of the bonds (here 0.8 × 0.8 cm², the flat part
of the face). The block has 540 interfaces and 2160 bonds.

**The law.** A bond is a spring (`knInnerBond`, `ktInnerBond`, per bond),
damped by `en2InnerBond`. The interface breaks, all its bonds at once, when
the elastic energy of its bonds exceeds `2 A Gc`, with `Gc = 50 J/m²`
(`gcInnerBond`): it is the energy needed to create the two faces of a crack.
After that, the two cells only interact by contact (`knContact 1 1`...), with
friction.

**What the block does.** `damage.txt` (`TrackDamage`) gives the broken
fraction of the interfaces: it rises to 0.567 within 0.7 milliseconds of
the impact, and no longer changes. `ClusterGranulo-conf20.txt` gives the
fraction, in number, of the fragments made of 1, 2, ..., 9 cells, and of 10
cells or more (the last class):

```
1 0.541667
2 0.208333
3 0.125
4 0
5 0
6 0
7 0.0416667
8 0
9 0.0138889
10 0.0694444
```

There are 72 fragments (0.0138889 = 1/72): 39 single cells, 15 pairs, 9
triplets, three fragments of 7 cells, one of 9, and five of 10 cells or more,
among them the top of the block.

**The energy balance.** `energies.sh` reads the cells and the interfaces of
`conf0` and of `conf20`:

```
kinetic energy before the impact        26.850 J
energy to break all the interfaces       3.456 J
energy spent in the cracks (conf20)      1.958 J  (7.3% of the kinetic energy)
kinetic energy left (conf20)             5.104 J
```

Only 7% of the kinetic energy opens cracks. The fragments carry away 5.1 J
(they are still flying at `t = 0.01 s`), and the rest is dissipated by the
damping of the bonds and of the contacts, and by friction.

### What to look at

- In `see`, the cone crushed under the corner and the fragments (the key `d`
  hides the plate). The keys `+` and `-` go to the next and the previous
  conf-files.
- The damage against time:

  ```sh
  gnuplot -p -e "plot 'damage.txt' using 1:2 with lines title 'damage'"
  ```
- The kinetic energy, in `kineticEnergy.txt` (`gnuplot checkplots.txt`).

### Exercises

1. Change the velocity in `setAllVelocities`: the damage is 0.06 at 2.5 m/s,
   0.17 at 5 m/s, 0.57 at 10 m/s and 0.93 at 20 m/s.
2. Make the material tougher: with `gcInnerBond 1 1 200`, the damage at
   10 m/s falls to 0.22.
3. Damp the bonds less (`en2InnerBond 1 1 0.9`): the elastic waves that run
   through the block are less damped, and break more interfaces (damage 0.70).
4. Make the cells **soft** (Rockable compiled with
   `ROCKABLE_ENABLE_SOFT_PARTICLES`, the default): add, before `stickBCM`,

   ```
   useSoftParticles 1e9 0.25
   useViscousTime 2e-5
   ```

   The cells now deform (Young's modulus 1 GPa, Poisson's ratio 0.25), and
   take a part of the energy: the damage is 0.47 (0.41 with `1e8`). Remove
   `useViscousTime`: the interfaces are now stiffer than the cells, the update
   of their deformation is unstable, and the computation diverges (the kinetic
   energy reaches 10²¹ J). A viscosity of 20 time steps
   (`2e-5 = 20 × dt`) is enough to stabilise it.
5. Drop the block on a **face**: in `genBlock.txt`, replace the two lines
   `transformation:...` by `transformation:add_translation 0 0.035 0`, run
   `generator genBlock.txt`, and put its lines in place of the cells of
   `input.txt`. The compression wave reflects on the top face as a tension,
   and the block splits into horizontal layers (*spalling*): the damage is
   only 0.13, and the top of the block flies back up.
6. Follow the fragmentation in time: `firstConf 0` and `stepConf 2` in
   `fragments.txt` give one file `ClusterGranulo-conf<n>.txt` per conf-file.
7. Scatter the toughness: add, after `stickBCM`,

   ```
   setVariableStickParams Gc inner 50 4 0
   ```

   Each interface now draws its own `Gc` from a Weibull law of scale 50 J/m²
   and modulus 4 (between 9 and 78 J/m² here, 45 J/m² on average). The weak
   interfaces break first, and the damage rises to 0.64; with a modulus of 2
   (a wider scatter), to 0.70. The last `0` makes the draw reproducible.

### Going further

- [`stickVerticesInClusters`](../stickVerticesInClusters): the other breakable
  law, `StickedLinks`, with a force threshold instead of an energy.
- [`softParticles`](../softParticles): soft cells, a bonded beam, a ball of
  soft cells that can break (`softBallImpact`).
- `setVariableStickParams` scatters the other parameters of the interfaces as
  well (stiffness, strength): [preProcessing](https://richefeu.github.io/rockable/preProcessing.html).
