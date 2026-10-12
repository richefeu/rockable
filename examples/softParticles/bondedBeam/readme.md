## Beam of bonded soft cubes

A simply supported beam made of 21 x 3 sphero-cubes (edge 1 cm, radius 1 mm),
glued face to face by `stickBCM` and bent by its own weight. It is a shorter
version (span 0.2 m instead of 0.6 m) of the beam of the DCEM paper, with the
same material: cells with `E = 1e8` and `nu = 0.3`, bonds with
`kn = kt = 1e5` per bond. It runs in about 5 seconds.

Rockable has to be compiled with `-DROCKABLE_ENABLE_SOFT_PARTICLES=ON`, and the
computation is run with

```sh
rockable input.txt
```

The left end of each row is fixed, and the right end can only slide along `x`
(`drivingSystem.txt`). The bonds are critically damped (`en2InnerBond 0 0 0.0`),
so the beam settles in about 0.2 s. The vertical position of the cube at
mid-span is recorded in `midspan.txt` (third column).

### What to check

| Cells                                          | Mid-span deflection |
|------------------------------------------------|---------------------|
| rigid (line `defineElasticParticles` commented) | 0.172 mm            |
| soft (as given)                                | 0.212 mm            |

The deformation of the cells makes the beam softer: the ratio of the effective
moduli is `0.172 / 0.212 = 0.81`. In the paper, this ratio is predicted from the
stiffness partition ratio

    eta = Kn (1 - 2 nu) / (E h) = 4e5 x 0.4 / (1e8 x 0.01) = 0.16

(`Kn = 4 kn` for the 4 bonds of a face, `h` the edge of a cell) by

    E_eff / E_eff_rigid = ((N^2 - 1) zeta (1 - 2 nu) / (eta + zeta (1 - 2 nu)) + chi) / (N^2 - 1 + chi)

with `N = 3` cells across the height, `zeta = 0.976` and `chi = 1.92`, which
gives 0.77. The difference (6 %) is of the order of the accuracy reported in the
paper for the longer beam (4 %).

The stiffnesses of the bonds are given for each bond. They could also be given
for a whole interface (`bondStiffnessPerInterface 1`, with `kn = 4e5`), with
the same result.
