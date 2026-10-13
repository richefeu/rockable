## 18. bunnyFun 4: jelly bunnies

<p align="center"><img src="animation.gif" alt="animation of the example"/></p>

### Goal

Three bunnies fall side by side from 35 cm onto a floor. They are **soft**:
each one deforms, and they differ only by their Young modulus, from a rubber
(1 MPa) to a soft jelly (10 kPa). The softest is squashed by 40% at the
impact. The soft particles of Rockable (the DCEM) are used here on a shape
that is neither a cube nor a cell: the mesh of the bunny, with its 157
vertices.

### Notions

- **Soft particles**: a homogeneous deformation of each particle, its
  **deformation gradient** `F` ([softParticles](https://richefeu.github.io/rockable/softParticles.html)).
- The keywords **`defineElasticParticles`** (a material for some particles)
  and **`useViscousTime`** (a viscosity that keeps the update stable).
- The **stability limit** of the soft particles: contacts stiffer than the
  particle need a viscosity.
- The section **`SoftParticles`** of the conf-files, read by a script.

Rockable must be compiled with `ROCKABLE_ENABLE_SOFT_PARTICLES` (the default):
`rockable -b` lists the compilation options.

### Files

| File            | Content                                                             |
|-----------------|---------------------------------------------------------------------|
| `../bunny.stl`  | The Stanford bunny, 309 triangles (see [`../readme.md`](../readme.md)) |
| `genShapes.txt` | Commands of `generator` that write `floor.shp`                      |
| `input.txt`     | The conf-file: the floor and the three soft bunnies                  |
| `squash.sh`     | The vertical stretch of each bunny, along the computation (sh + awk) |

### Run it

```sh
stl2shape -r 0.002 -c -n Bunny -o bunny.shp ../bunny.stl
shapeSurvey bunny.shp -P bunny.shp
generator genShapes.txt                      # writes floor.shp
cat floor.shp bunny.shp > shapes.shp
rockable input.txt                           # less than a second, conf0 to conf200
see conf100
sh squash.sh > squash.txt
```

### What happens

**The soft bunnies.** In `input.txt`, after the line `Particles`, three lines
make the bunnies (ranks 1, 2 and 3) elastic:

```
defineElasticParticles 1e6 0.3 1 1
defineElasticParticles 3e4 0.3 1 2
defineElasticParticles 1e4 0.3 1 3
useViscousTime 5e-3
```

the Young modulus, the Poisson ratio, the number of particles and their
ranks. A soft particle carries a deformation gradient `F`: each vertex `v` of
the shape is placed at `x + Q (F h v)`. The deformation is the same over the
whole bunny (it is *homogeneous*): a bunny can be squashed, stretched or
sheared, but not bent. At each step, `F` follows the stress that the contacts
apply to the particle, through the elasticity of the material.

**The stability.** A soft particle is updated explicitly, and the update is
stable only when the contacts are softer than the particle, about
`kn < E h`. Here `kn = 3e4 N/m` per contact, and a bunny of 10 cm with
`E = 1e4 Pa` gives `E h = 1e3 N/m`: the contacts are far too stiff. The line
`useViscousTime 5e-3` gives the bunnies a viscosity, over 5 ms (50 time
steps), which makes the update stable without changing the response at rest
(exercise 1 removes most of it).

**The impact.** The bunnies fall for 0.25 s, and hit the floor at 2.5 m/s.
`squash.sh` turns `F` into the frame of the floor and gives its vertical
component `F_yy`: 1 at rest, 0.6 for a bunny squashed by 40%.

```
sh squash.sh | tail -1
# smallest: 0.9874 0.7926 0.5817
```

The rubber bunny is squashed by 1.3%, the medium one by 21%, the jelly by 42%.
They then recover their shape, bounce a little, roll over, and rest on the
floor; the two softest stay squashed by about 1% by their own weight.

**What the model does not do.** `F_yy` never goes above 1: the bunnies are
not stretched when they bounce back, and they do not wobble like a jelly. In
this model, the deformation has no inertia of its own: it follows the forces
of the contacts, with the delay of the viscosity. A wobble would need the
mass of the particle to be spread over its vertices.

### What to look at

- In `see`, the bunnies are drawn deformed; `conf102` to `conf112` show the
  impact. The keys `4` and `5` colour the soft particles by their change of
  volume `det F - 1` and by their mean stress. The skins of `shape2mesh` cannot
  follow the deformation: soft particles are always drawn with their spheres,
  tubes and slabs.
- The squash against time:

  ```sh
  sh squash.sh > squash.txt
  gnuplot -p -e "plot for [c=2:4] 'squash.txt' using 1:c with lines title sprintf('bunny %d', c-1)"
  ```
- In a conf-file, the section `SoftParticles`: the material, `F`, the stress
  and the strain of each soft particle, in the frame of the particle.

### Exercises

1. **The stability limit.** Write `useViscousTime 1e-3`: the jelly bunny
   explodes at the impact (`F_yy` reaches -507, the kinetic energy becomes
   `nan`), and the medium one, less damped, is squashed by 32%. The
   viscosity has to be a few tens of time steps.
2. **A more viscous jelly.** `useViscousTime 1.5e-2`: the bunnies resist the
   impact more, and are squashed by 11% and 26% only.
3. **The Poisson ratio.** Write `0.45` instead of `0.3` in the three lines
   `defineElasticParticles`. The squash hardly changes (40% for the jelly),
   but the volume does: `det F - 1` reaches -20% with `0.3`, -5% with `0.45`
   (key `4` in `see`, or the section `SoftParticles`). A nearly incompressible
   bunny widens instead of shrinking.
4. **Rigid bunnies.** Remove the four lines that make the bunnies soft: the
   three bunnies now do exactly the same, and bounce higher (their centre
   rises to 11 cm, against 9.2 cm for the rubber bunny, 4.3 and 5.8 cm for the
   soft ones): the viscosity of the soft bunnies dissipates energy at the
   impact.

### Going further

- [`softParticles`](../../softParticles): soft cubes, a bonded beam, a ball of
  soft cells that breaks.
- The DCEM model, its stability and its keywords:
  [softParticles](https://richefeu.github.io/rockable/softParticles.html).
- [`bunnyShot`](../bunnyShot): the bunny cut into bonded cells.
