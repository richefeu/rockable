## Viscoelastic cube between two platens

The same test as `elasticCube`, but with contacts a hundred times stiffer
(`knContact 0 1 1e7`). With the elastic update, the cube then oscillates from
one time step to the next and never settles (see the readme of `elasticCube`).

The cube is made viscoelastic (Kelvin-Voigt) with

```
useSoftParticles 1e6 0.0
useViscousTime 1e-2
```

that is, a viscosity proportional to the stiffness, `D = tau C`, with
`tau = 1e-2 = 100 dt`. The strain no longer follows the stress without delay:
it relaxes towards its elastic value over a time `tau`, which removes the
instability. The viscous stress vanishes at equilibrium, so the final state is
the elastic one: `F_yy = exp(sigma_yy / E)`, close to 0.99.

Rockable has to be compiled with `-DROCKABLE_ENABLE_SOFT_PARTICLES=ON`, and the
computation (less than a second) is run with

```sh
rockable input.txt
```

### Choosing tau

In the DCEM paper, the explicit update is shown to be stable when
`tau > (1 + xi) dt / 2`, where `xi` is the loop gain between the stiffness of the
interfaces and that of the particle. A viscosity of a few tens of time steps is
usually enough; a larger `tau` slows down the deformation of the particle.
