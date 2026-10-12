## Elastic cube between two platens

A soft cube (Young modulus 1e6, Poisson ratio 0) is squeezed between a fixed
bottom platen and a top platen pushed down by a force of 1e4. It runs in less
than a second.

Rockable has to be compiled with

```sh
cmake .. -DROCKABLE_ENABLE_SOFT_PARTICLES=ON
```

and the computation is run with

```sh
rockable input.txt
```

### What to check

The cube is the only soft particle (`useSoftParticles` makes soft all the
particles that are not driven). At equilibrium, its mean stress is

    sigma_yy = -F L / V

where `F = 1e4` is the force, `L` the distance between the contact points of
the two platens and `V` the deformed volume of the cube. With a zero Poisson
ratio, the cube only shortens along `y`, and the component `yy` of the
deformation gradient is `exp(sigma_yy / E)`, close to 0.99.

Both can be read in the section `SoftParticles` of the last conf-file:

```
SoftParticles 1
#id young poisson F(9) stress(9) strain(9)
2 1e+06 0 1 0 0 0 0.98999... 0 0 0 1 ...
```

### Stiffness of the contacts

The contact stiffness (`knContact 0 1 1e5`) is lower than the stiffness of the
body, which is about `E A / L = 1e6`. The deformation of a soft particle follows
its stress without delay, so a contact much stiffer than the body makes the
computation oscillate from one time step to the next: the body shortens so much
that it loses the contact, recovers its shape, touches again, and so on. With
`knContact 0 1 1e7`, this example no longer converges. A small viscosity of the
particle removes this problem: see the example `viscousCube`.
