## Impact of a ball of soft cells

A ball of 123 sphero-cubic cells (edge 0.3 mm), glued face to face by
`stickBCM`, hits a rigid plate at 5 m/s. The cells are very soft
(`E = 1 MPa`, `nu = 0.3`) and viscoelastic, so the ball flattens on the plate,
recovers its shape and rebounds. It runs in about 3 seconds.

This is a small version of the illustration of the DCEM paper (Figure "soft
cells"), with the same parameters but cubic cells instead of Voronoi cells. The
full case of the paper is in `../softBallImpactVoronoi` (about 12 minutes).

Rockable has to be compiled with `-DROCKABLE_ENABLE_SOFT_PARTICLES=ON`, and the
computation is run with

```sh
rockable input.txt
```

The units are mm, g and ms, so that stresses are in MPa and velocities in m/s.

### Why a viscosity

The interfaces are about as stiff as the cells, far beyond the stability limit
of the elastic update of the strain (see `../elasticCube`). The cells are
therefore given a viscosity proportional to their stiffness,

```
defineElasticParticles 1 0.3 123 1 2 ... 123
useViscousTime 1e-2
```

that is `D = tau C` with `tau = 100 dt`.

### What to check

- `lowestCell.txt` (`TrackSoftParticle`) follows the cell at the bottom of the
  ball: its volume ratio `det(F)` (second column) goes down to about 0.49 at
  `t = 0.15 ms`, and comes back to about 0.97 after the rebound. The small
  residual strain is held by the tangential bond forces, which are computed
  incrementally (the paper reports the same effect).
- The cells leave the plate at about 2.8 m/s: the restitution coefficient is
  about 0.55.
- In `see`, the keys `4` (volume change `det(F) - 1`) and `5` (mean stress) color
  the cells.

### Breakable interfaces

With `gcInnerBond 0 0 3e-3` instead of `1000`, about 40 % of the interfaces
break (`damage.txt`, `TrackDamage`): the lower cells are crushed against the
plate instead of rebounding.
