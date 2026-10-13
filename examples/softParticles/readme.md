## Soft particles (DCEM)

Examples of particles with a homogeneous deformation (the discrete
compliant-element method, DCEM). Rockable has to be compiled with

```sh
cmake .. -DROCKABLE_ENABLE_SOFT_PARTICLES=ON
```

| Folder                  | What it shows                                                         | Time   |
|-------------------------|-----------------------------------------------------------------------|--------|
| `elasticCube`           | An elastic cube squeezed between two platens; its stability limit     | < 1 s  |
| `viscousCube`           | The same with much stiffer contacts, stabilised by a viscosity        | < 1 s  |
| `bondedBeam`            | A beam of bonded soft cubes bent by its weight, compared with rigid cells | 5 s |
| `softBallImpact`        | A ball of 123 soft, viscous, bonded cells hitting a plate             | 3 s    |
| `softBallImpactVoronoi` | The same with the 200 Voronoi cells of the DCEM paper                 | 12 min (4 threads) |

Each folder has a `readme.md` that tells what to check in the results.

The tutorial [`fragmentation`](../fragmentation) breaks a block of bonded
cells (`BCM`) on impact, with rigid or soft cells.
