## 6. Breakable bonds: crushing a ring of glued sectors

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

A ring made of eight sectors, glued together, is crushed between two plates
until it breaks. The glue is made of breakable bonds between the vertices of
neighbouring sectors. This example was set as a benchmark for the PhD thesis
of Marta Stasiak.

### Notions

- The **cluster** number of a particle (second number after the shape name).
- The pre-processing command **`stickVerticesInClusters`**, which glues the
  close vertices of the particles of a same cluster
  ([preProcessing](https://richefeu.github.io/rockable/preProcessing.html)).
- The force law **`StickedLinks`** and the parameters of its bonds
  (`...InnerBond`) ([forceLaws](https://richefeu.github.io/rockable/forceLaws.html)).
- The **interfaces** saved in the conf-files.
- The data extractor **`TrackBody`**: the force on a plate.

### Files

| File                 | Content                                                          |
|----------------------|------------------------------------------------------------------|
| `input.txt`          | The conf-file: two plates and eight sectors                       |
| `sector.shp`         | The shapes `Plan` and `ShellSector8_N` (one eighth of a ring)     |
| `drivingSystem.txt`  | The top plate comes down at 1 cm/s                               |
| `dataExtractors.txt` | Records the position and the force of the two plates             |

### Run it

```sh
rockable input.txt           # about 2 seconds
see
```

### What happens

**The system.** The first two particles are the plates (`nDriven 2`), in
group 1. The eight sectors, in group 0, form a ring of about 17 mm standing
between them. All the sectors have the cluster number `0`:

```
ShellSector8_N 0 0 1 ...
               | |
               | cluster
               group
```

**The glue.** The last line of `input.txt`,

```
stickVerticesInClusters 0.0004
```

is run when the file is read: each pair of vertices of two particles of the
same cluster closer than 0.4 mm is glued by a bond. Here, each sector shares
four vertices with each of its two neighbours: the ring holds 8 interfaces of
4 bonds. The glue is put between the vertices only, and its normal direction is
the branch between the centres of the two particles (the normal of the faces
is not used). The relative displacements are set to zero when a bond is made,
so that all the forces start at zero.

**The law.** `forceLaw StickedLinks` handles two kinds of interactions: the
usual frictional contacts (`...Contact` parameters, for the pairs of groups
`0 1` and `0 0`), and the bonds between the vertices of a same cluster
(`...InnerBond` parameters, for the pair of groups `0 0`):

```
knInnerBond 0 0 0.8e7      normal stiffness of a bond
en2InnerBond 0 0 0.1       restitution (damping) of a bond
ktInnerBond 0 0 0.8e7      tangential stiffness
fn0InnerBond 0 0 85.       normal force at failure (tension)
ft0InnerBond 0 0 170.0     tangential force at failure
powInnerBond 0 0 2.0       exponent of the failure criterion
```

A bond fails when

```
(|ft| / ft0)^pow + (|mom| / mom0)^pow - fn / fn0 > 1
```

where `fn` is positive in compression and negative in tension, and `mom` the
moment carried by the bond (zero here: there is no `krInnerBond`). A
compression therefore strengthens the bond. When one bond fails, **all the
bonds of its interface** are broken, and the two sectors then only interact by
frictional contacts.

**The loading.** `drivingSystem.txt` holds `Control _y_Vel_ 1 -0.01`: the top
plate comes down at 1 cm/s. The bottom plate is fixed. There is no gravity, and
no damping.

The force on the top plate grows linearly up to about **121 N**, at a
displacement of 0.1 mm (`t = 0.0095 s`): two interfaces break, and the force
drops. It rises again (about 85 N at 0.2 mm), then two more interfaces break
around `t = 0.023 s`. The ring is then in four pieces of two sectors, which are
thrown sideways: the force falls to zero.

### What to look at

- In `see`, the sectors separate at the interfaces that break.
- The force on the top plate (column 16 of `plan_up.txt`) against the
  displacement of the plate (column 3), in millimetres:

  ```sh
  gnuplot -p -e "plot 'plan_up.txt' using ((0.0174844-\$3)*1000):16 with lines title 'force (N)'"
  ```

  (`0.0174844` is the initial position of the plate.) `extractedDataDoc.txt`
  lists the columns of the files written by the data extractors.
- The number of interfaces still glued is saved in each conf-file:

  ```sh
  for i in $(seq 0 80); do echo "conf$i $(grep '^Interfaces' conf$i)"; done
  ```

  It goes from 8 to 6 between `conf18` and `conf20`, then to 4.

### Exercises

1. Double the tensile strength of the bonds (`fn0InnerBond 0 0 170.`): the
   peak force doubles too (about 241 N at 0.2 mm). The ring breaks in tension.
2. Set `powInnerBond 0 0 1.0`. Nothing changes: the exponent only applies to
   the tangential force and to the moment, and the bonds fail in tension.
3. Give a different cluster number to each sector (`0` to `7`): no bond is
   made (the conf-files hold `Interfaces 0`), and the ring opens at once,
   without any resistance.
4. Crush ten times faster (`Control _y_Vel_ 1 -0.1`): the peak is a little
   lower (110 N), and four interfaces break at once.
5. Set `tmax 0.1`, and follow the pieces with `see` once the ring is broken.

### Going further

- [`DominoFun`](../DominoFun): a wall of glued bricks hit by dominoes.
- [`softParticles`](../softParticles): the `BCM` law, which glues the faces of
  deformable particles.
- The bonds and their failure criterion:
  [forceLaws](https://richefeu.github.io/rockable/forceLaws.html).
