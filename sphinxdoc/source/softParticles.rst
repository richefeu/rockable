.. _softParticles:

Deformable particles (DCEM)
===========================

By default a Rockable particle is rigid: it has six degrees of freedom, three of
translation and three of rotation. The *soft particles* extension gives some of
the particles a **uniform deformation**, so that they can change shape under
load without being discretised into a mesh. This is the discrete
compliant-element method (DCEM), described in the paper *The Discrete
Compliant-Element Method in Rockable framework*. It was developed in the
postdoctoral work of *Mukesh Singh Bisht*.

The intent is not to resolve a full kinematic field with many degrees of
freedom, but to give each particle the few that matter when the deformation of
the particles is part of the physics: soft grains, cells of a bonded material
whose fracture depends on the elastic energy they store, polymers or biological
tissues.

.. contents::
   :local:
   :depth: 2


Compilation
-----------

The feature is enabled by the CMake option ``ROCKABLE_ENABLE_SOFT_PARTICLES``,
which is ``ON`` by default. It can be disabled when it is not needed:

.. code-block:: sh

   cmake .. -DROCKABLE_ENABLE_SOFT_PARTICLES=OFF

Without it, the keywords below are accepted but only print
``SOFT_PARTICLES_NOT_ENABLED when Rockable was compiled``.

.. tip::

   ``rockable -b`` prints ``ROCKABLE_ENABLE_SOFT_PARTICLES`` among the
   compilation options, which is the quickest way to tell whether the binary at
   hand supports the feature.


The model
---------

Kinematics
^^^^^^^^^^

Each soft particle :math:`i` carries a deformation gradient
:math:`\mathbf{F}_i`, expressed in the frame of the body. A vertex
:math:`\underline{v}` of the shape is placed at

.. math::

   \underline{x} = \underline{x}_i + \mathbf{Q}_i \left( \mathbf{F}_i \, h_i \underline{v} \right)

where :math:`\underline{x}_i` is the position of the particle,
:math:`\mathbf{Q}_i` its orientation and :math:`h_i` its homothety. The
deformation acts on the polyhedral core only: the Minkowski radius (the
thickness of the sphero-polyhedron) is **not** deformed. The contact detection,
the oriented bounding boxes, ``see`` and ``conftovtk`` all use the deformed
vertices. The mass and the inertia of the particle are not changed by the
deformation.

Stress
^^^^^^

At each computation of the forces, the mean (Love--Weber) stress in the
particle is computed from the forces :math:`\underline{f}_c` it receives at its
contacts and bonds :math:`c`, at the points :math:`\underline{x}_c`:

.. math::

   \boldsymbol{\sigma}_i = \frac{1}{V_i} \sum_c \underline{f}_c \otimes \left( \underline{x}_c - \underline{x}_i \right),
   \qquad V_i = h_i^3 \, V^\text{shape} \det \mathbf{F}_i

symmetrised, and expressed in the frame of the body. With this convention, a
compressive stress is negative.

Elastic material
^^^^^^^^^^^^^^^^

The material is linear, isotropic and elastic, with a Young modulus :math:`E`
and a Poisson ratio :math:`\nu`, and :math:`\mathbb{C}^{-1}` its compliance. The
increment of stress since the last update gives an increment of strain, which is
applied exactly to the deformation gradient:

.. math::

   \Delta \boldsymbol{\varepsilon} = \mathbb{C}^{-1} : \left( \boldsymbol{\sigma}^{n} - \boldsymbol{\sigma}^{n-1} \right),
   \qquad
   \mathbf{F}^{n} = \exp \left( \Delta \boldsymbol{\varepsilon} \right) \mathbf{F}^{n-1}

The strain follows the stress without delay, so the update does not depend on
the time step.

Viscoelastic material
^^^^^^^^^^^^^^^^^^^^^

A Kelvin--Voigt material adds a viscosity :math:`\mathbb{D}`,
:math:`\boldsymbol{\sigma} = \mathbb{C} : \boldsymbol{\varepsilon} + \mathbb{D} : \dot{\boldsymbol{\varepsilon}}`.
The strain is then integrated explicitly,

.. math::

   \boldsymbol{\varepsilon}^{n} = \boldsymbol{\varepsilon}^{n-1}
   + \Delta t \, \mathbb{D}^{-1} : \left( \boldsymbol{\sigma}^{n} - \mathbb{C} : \boldsymbol{\varepsilon}^{n-1} \right),
   \qquad
   \mathbf{F}^{n} = \exp \left( \boldsymbol{\varepsilon}^{n} - \boldsymbol{\varepsilon}^{n-1} \right) \mathbf{F}^{n-1}

:math:`\mathbb{D}^{-1}` has the same form as the compliance, built from a shear
viscosity :math:`\eta_s` and a bulk viscosity :math:`\eta_b` (in place of
:math:`G` and :math:`K`). The simplest choice is a viscosity proportional to the
stiffness, :math:`\mathbb{D} = \tau \, \mathbb{C}`, that is
:math:`\eta_s = \tau G` and :math:`\eta_b = \tau K`: the strain then relaxes
towards its elastic value over the time :math:`\tau`. The viscous stress
vanishes at equilibrium, so the quasi-static response is the elastic one.

.. warning::

   The viscoelastic strain is integrated over :math:`\Delta t` each time the
   forces are computed. The integrators ``Beeman`` and ``RungeKutta4`` compute
   them more than once per time step, so the viscoelastic particles are only
   correct with ``Euler`` and ``velocityVerlet`` (the default). A warning is
   printed otherwise.


Stability of the elastic update
-------------------------------

The deformation gradient is updated once per computation of the forces,
without any relaxation. This explicit iteration converges only if the
interfaces of a particle (contacts or bonds) are **softer** than the particle
itself. Otherwise, the particle shortens so much within one step that it loses
its contacts, recovers its shape, touches again, and so on: the computation
oscillates from one time step to the next, whatever the time step. In the DCEM
paper, the limit is written with the stiffness partition ratio

.. math::

   \eta = \frac{K_n (1 - 2\nu)}{E\, h} < \eta_\text{lim} = (1 - 2\nu) \frac{\zeta h}{\ell}

where :math:`K_n` is the normal stiffness of a whole interface (the sum of its
bonds or contacts), :math:`h` the size of the particle, :math:`\ell` the size of
its core and :math:`\zeta` the ratio of its volume to :math:`h^3`. In short, the
contact stiffness must stay below about :math:`E h`.

A small viscosity of the particles removes this limit, without changing the
quasi-static response. The paper finds the update stable when

.. math::

   \tau > \frac{1 + \xi}{2} \Delta t

where :math:`\xi` is the loop gain between the stiffness of the interfaces and
that of the particle. A viscosity of a few tens of time steps is usually enough,
and is recommended whenever interfaces can break: a broken bond turns into a
contact, which can be much stiffer.


Keywords
--------

The particles have to be read before these keywords, which therefore come after
``Particles`` in the input file.

- ``useSoftParticles`` (*double*) **Young** (*double*) **Poisson**

  All the free particles (the ones that are not driven) become elastic, with the
  given material.

- ``defineElasticParticles`` (*double*) **Young** (*double*) **Poisson** (*int*) **number** (*int*) **id1** **id2** ...

  The listed particles (identified by their rank in ``Particles``) become
  elastic, with the given material. A particle that is already soft keeps its
  deformation and only changes its material.

- ``useViscousTime`` (*double*) **tau**

  All the soft particles become viscoelastic, with
  :math:`\mathbb{D} = \tau \mathbb{C}`. It must come after the elasticity has
  been set.

- ``defineViscousTime`` (*double*) **tau** (*int*) **number** (*int*) **id1** **id2** ...

  The same, for the listed particles.

- ``defineViscousParameters`` (*double*) **shearViscosity** (*double*) **bulkViscosity** (*int*) **number** (*int*) **id1** **id2** ...

  The listed soft particles become viscoelastic, with the given shear and bulk
  viscosities. This is the keyword used in the simulations of the DCEM paper,
  with :math:`\eta_s = \tau G` and :math:`\eta_b = \tau K`.

.. code-block:: text
   :caption: input.txt

   nDriven 2
   Particles 3
   Platen_b 1 1 1   0 -1.1 0   0 0 0   0 0 0   1 0 0 0   0 0 0   0 0 0
   Platen_t 1 1 1   0  0.1 0   0 0 0   0 0 0   1 0 0 0   0 0 0   0 0 0
   c        0 0 1   0 -0.5 0   0 0 0   0 0 0   1 0 0 0   0 0 0   0 0 0

   useSoftParticles 1e6 0.3
   useViscousTime 1e-2

.. note::

   The stiffness of the **contacts** is a different matter, set by
   ``knContact`` and ``ktContact`` as usual: a soft particle in a stiff contact
   and a stiff particle in a soft contact are two different models, and both are
   available within the stability limit above.


Bonded soft cells
-----------------

Soft particles can be glued together by the pre-processing command
``stickBCM`` and the force law ``BCM`` (see :ref:`Force-laws`), to model a
material made of cells that deform and whose interfaces can break. In the DCEM
paper, the stiffnesses of the bonds are given for **each bond**, which is the
default; ``bondStiffnessPerInterface 1`` gives them for a whole interface
instead.

The bonds and the cells represent the same material: their stiffnesses must be
chosen together, and both count in the effective modulus of the assembly. The
paper gives the relations to use.


Conf-files and restart
----------------------

The state of the soft particles is written in the conf-files, right after the
particles:

.. code-block:: text

   SoftParticles 1
   #id young poisson shearViscosity bulkViscosity F(9) stress(9) strain(9)
   2 1e+06 0 0 0 1 0 0 0 0.98999 0 0 0 1 ...

The material, the deformation gradient, and the stress and strain (in the frame
of the body) are kept, so a computation restarted from a conf-file follows the
uninterrupted one. The viscosities are zero for an elastic particle. The
section ``Transformations``, with only the deformation gradients, written by the
first versions of the feature, is still read.


Output and visualisation
------------------------

- The data extractor ``TrackSoftParticle`` records the deformation gradient,
  the stress and the strain of one particle, and ``TrackStress`` the component
  :math:`\sigma_{xx}` of many (see :ref:`dataExtractors`).
- In ``see``, the soft particles are drawn deformed, and the keys ``4`` and
  ``5`` color them by their volume change :math:`\det \mathbf{F} - 1` and by
  their mean stress :math:`\operatorname{tr} \boldsymbol{\sigma} / 3`.
- ``conftovtk`` writes the deformed shapes.


Examples
--------

The folder ``examples/softParticles`` holds examples that run in a few seconds,
each with a ``readme.md`` that tells what to check:

- ``elasticCube``: an elastic cube squeezed between two platens; the computed
  strain is :math:`\sigma / E`.
- ``viscousCube``: the same with contacts a hundred times stiffer, beyond the
  stability limit of the elastic update, stabilised by a viscosity.
- ``bondedBeam``: a beam of bonded soft cubes bent by its weight, compared with
  rigid cells.
- ``softBallImpact``: a ball of soft, viscous, bonded cells hitting a plate,
  with unbreakable or breakable interfaces.
- ``softBallImpactVoronoi``: the illustration of the DCEM paper, with 200
  Voronoi cells (about 12 minutes).
