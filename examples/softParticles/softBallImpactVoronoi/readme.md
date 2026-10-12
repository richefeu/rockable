## Impact of a ball of soft Voronoi cells (DCEM paper)

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

The illustration of the DCEM paper (Figure "soft cells", case (a)): a sphere
made of 200 Voronoi cells, glued by `stickBCM`, hits a rigid plate at 5 m/s.
The cells are very soft (`E = 1 MPa`, `nu = 0.3`) and viscoelastic
(`tau = 100 dt`, given here by `defineViscousParameters` with
`shearViscosity = tau G` and `bulkViscosity = tau K`).

This is a **long** example: about 12 minutes with 4 threads,

```sh
rockable -j 4 input.txt
```

with Rockable compiled with `-DROCKABLE_ENABLE_SOFT_PARTICLES=ON`. The units
are mm, g and ms. A quicker version, with cubic cells, is in
`../softBallImpact`.

### What to check

The results of the paper are recovered:

- the cells near the contact lose up to about 40 % of their volume (the
  smallest `det(F)` is about 0.62 at `t = 0.1 ms`);
- the ball rebounds at 4.34 m/s, that is a restitution coefficient of 0.87.

With `gcInnerBond 0 0 3e-2` (3 J/m^2 in the units of the paper), the interfaces
near the contact break, and the ball ends as a flattened cap on its own debris
(damage `D = 0.46` in the paper).
