.. _Force-laws:

Force laws
==========

We provide here the interaction models existing in ``Rockable``.
The equations are explained with the associated parameters.
The parameters ``[PARAMETER]`` are to be defined in the input conf-file in the following way: 

- ``[PARAMETER]`` (*int*) **firstGroupNumber** (*int*) **secondGroupNumber** (*double*) **value**

For example, to set the value :math:`10^8` to the normal stiffness :math:`k_n` between the elements that belong to the group number ``0`` and the elements that belong to the group number ``2``, we will write:

.. code-block:: text
   :caption: input.txt
	 
	 ...
   knContact 0 2 1e8
	 ...


Default model (keyword ``Default``)
-----------------------------------

This is the force law to be used in most cases. In short, the model includes a
linear elastic normal force, a normal viscous damping, a Coulomb friction, and
a rolling resistance moment. There is no cohesion: the interaction exists only
while the bodies overlap (:math:`d_n \leq 0`), and all its forces are reset to
zero as soon as the contact is lost.

1. **Normal component**

  a. The elastic part of the normal contact force is

  .. math::

     f_n^\text{elas} = -k_n d_n

  where :math:`d_n \leq 0` is the normal distance, so :math:`f_n^\text{elas} \geq 0`. The value of :math:`k_n` is set with the keyword ``knContact``.

  b. The viscous part of the normal contact force is

  .. math::

     f_n^\text{visc} = 2 \alpha_n \sqrt{k_n m_\text{eff}} \, v_n

  where :math:`v_n` is the relative normal velocity, :math:`m_\text{eff}=(m_i m_j)/(m_i+m_j)` is the effective mass, and :math:`\alpha_n \in [0, 1[` is the rate of normal viscous damping. There are two ways to set the value of :math:`\alpha_n`:

    i. The keyword ``en2Contact``, which is the squared normal restitution coefficient (:math:`e_n^2`). The viscous damping rate is then

    .. math::

       \alpha_n = \frac{- \ln e_n}{\sqrt{\ln^2 e_n + \pi^2}}

    ii. The keyword ``en2ContactFromViscRate``, which gives :math:`\alpha_n` directly; the corresponding restitution is

    .. math:: 

       e_n^2 = \exp \left(-\frac{2 \alpha_n \pi}{\sqrt{1 - \alpha_n^2}}\right)

  .. note:: 

     Because of the viscous part, the total normal force :math:`f_n = f_n^\text{elas} + f_n^\text{visc}` could be negative. It is restricted to remain positive or zero.

  .. note::

     With ``preventCrossingLength`` :math:`L > 0`, the elastic part is
     multiplied by :math:`L / (L + d_n)`, which stiffens the repulsion at large
     overlap and prevents the bodies from crossing each other. The factor is
     capped at 100. This holds for ``Default`` and ``GeoVisc``.


2. **Tangential component**

The tangential force is accumulated incrementally from the relative tangential
velocity :math:`\underline{v}_t`, with the tangential stiffness :math:`k_t`
(keyword ``ktContact``), and projected onto the Coulomb cone defined by the
friction coefficient :math:`\mu` (keyword ``muContact``) and the **elastic**
part of the normal force:

.. math::

   \underline{f}_t \leftarrow \underline{f}_t + k_t \underline{v}_t \Delta t,
   \qquad
   \left\Vert \underline{f}_t \right\Vert \leq \mu f_n^\text{elas}

Compiling with ``ALLOW_NEGATIVE_THRESHOLDS`` uses the total normal force
:math:`f_n` instead of :math:`f_n^\text{elas}` in the threshold.


3. **Moment vector**

When the rolling stiffness :math:`k_r` (keyword ``krContact``) is strictly
positive, a resistant moment is accumulated from the relative rotation velocity,
and capped by the rolling resistance :math:`\mu_r` (keyword ``murContact``),
which is a **length**, times the normal force:

.. math::

   \underline{M} \leftarrow \underline{M} + k_r \left( \underline{\omega}_j - \underline{\omega}_i \right) \Delta t,
   \qquad
   \left\Vert \underline{M} \right\Vert \leq \mu_r f_n

.. note:: When :math:`k_r = 0`, no moment is computed.

.. note::

   When a contact partnership is active (keyword ``ContactPartnership``), the
   stiffnesses :math:`k_n`, :math:`k_t`, :math:`k_r` of each sub-contact are
   multiplied by its weight :math:`w`, and the damping coefficient by
   :math:`\sqrt{w}`.


Law for rock avalanches (keyword ``Avalanche``)
-----------------------------------------------

This is historically the first law that has been implemented in ``Rockable``
(actually in ``DEMbox``, its ancestor), to simulate rock avalanches at
Laboratoire 3SR. It differs from ``Default`` by the way the energy is
dissipated at the contacts: there is no viscous damping, but a hysteretic
normal force.

1. **Normal component**

   The normal force is computed incrementally while the bodies get closer
   (loading), and follows a softer unloading branch when they separate:

   .. math::

      f_n \leftarrow
      \begin{cases}
        f_n - k_n \, \Delta d_n & \text{if } \Delta d_n < 0 \text{ (loading)} \\
        - e_n^2 \, k_n \, d_n & \text{if } \Delta d_n > 0 \text{ (unloading)}
      \end{cases}
      \qquad f_n \geq 0

   where :math:`\Delta d_n` is the change of the normal distance since the
   previous step. The unloading stiffness is :math:`e_n^2 k_n`, so that the
   ratio of the energy given back to the energy stored is :math:`e_n^2`
   (keyword ``en2Contact``). The force is kept positive or zero.

2. **Tangential component**

   As in ``Default``, accumulated incrementally with :math:`k_t` and capped by
   the Coulomb cone, here with the total normal force:
   :math:`\Vert \underline{f}_t \Vert \leq \mu f_n`.

3. **Resistant moment**

   Accumulated as in ``Default``,
   :math:`\underline{M} \leftarrow \underline{M} + k_r (\underline{\omega}_j - \underline{\omega}_i) \Delta t`
   (even when :math:`k_r = 0`, which then gives no moment), and capped by

   .. math::

      \left\Vert \underline{M} \right\Vert \leq \mu_r \, \ell \, f_n

   where :math:`\mu_r` (keyword ``murContact``) is here **dimensionless**, and
   :math:`\ell` is the distance from the contact point to the rotation axis of
   the free body (a rock block), that is, the length of the branch vector from
   its centre to the contact point once its component along the rotation
   vector is removed. When body :math:`i` is driven, the free body is
   :math:`j`; otherwise it is :math:`i`. When the free body does not rotate at
   all, the whole branch vector is used.

.. note::

   The contact partnership weights :math:`k_n`, :math:`k_t` and :math:`k_r` as
   in ``Default``. ``FT_CORR`` applies the objectivity correction to the
   tangential force, as in the other laws.



Breakable elastic solid bonds (keyword ``StickedLinks``)
--------------------------------------------------------

This law has been introduced since the PhD work of *Marta Stasiak*.

Two states can exist: cohesive bond and contact.

1. **Cohesive bond**

    The parameters are those of the interface when ``ParamsInInterfaces`` is 1.
    Otherwise, they are read in the group-pair table, as ``*InnerBond`` when
    the two bodies share the same cluster and ``*OuterBond`` otherwise.

    a. **Normal component**

       The normal force (:math:`f_n = f_n^\text{elas} + f_n^\text{visc}`) between the particles is calculated using
       Hooke's elastic contact, combined with viscous damping:

       .. math::

          f_n^\text{elas} = -k_n (d_n - d_n^0)
        
       where :math:`d_n^0` is the normal distance of the interface when it was
       created. The parameter :math:`k_n` is set with the keyword
       ``knInnerBond`` or ``knOuterBond``. The viscous normal force is:

       .. math::

          f_n^\text{visc} = 2 \alpha_n \sqrt{k_n m_\text{eff}} \, v_n	

       where :math:`\alpha_n` is indirectly set with the keyword ``en2InnerBond`` or ``en2OuterBond``.
       The normal force is not clamped: a bond can carry a tensile force.

    b. **Tangential component**

       The tangential force is the sum of an elastic part, incrementally
       calculated from the relative tangential velocity :math:`\underline{v}_t`
       and the tangential stiffness :math:`k_t` (keyword ``ktInnerBond`` or
       ``ktOuterBond``), and of a viscous part, with the damping coefficient of
       the normal force:

       .. math::

           \underline{f}_t^e \leftarrow \underline{f}_t^e + k_t \underline{v}_t \Delta t,
           \qquad
           \underline{f}_t = \underline{f}_t^e + 2 \alpha_n \sqrt{k_n m_\text{eff}} \, \underline{v}_t

       As in ``BCM``, the viscous part of the previous step is removed before
       the elastic increment is added, so that it is not cumulated.

       .. warning::

          The viscous part is explicit, and :math:`m_\text{eff}` is built from
          the masses of the bodies, whereas a bond placed at a vertex also
          acts on their **rotation**. For flat or slender bodies held by many
          bonds, the effective inertia seen by a bond is much smaller than
          :math:`m_\text{eff}`, and a strong damping (a small ``en2*Bond``,
          close to 0) makes the integration unstable: the angular velocities
          grow exponentially from round-off, and the bonds break within a few
          hundred time steps, without any load. A larger ``en2*Bond`` or a
          smaller time step removes the problem. In ``examples/DominoFun``
          (bricks of 5 x 40 x 20 mm, ``kn = kt = 1e7``, ``dt = 1e-6``),
          ``en2InnerBond 0.02`` breaks the wall at once, ``0.5`` does not.

          Before the viscous part stopped being cumulated, this instability
          did not show: the cumulated term acted as an extra stiffness, not
          as a damping.

    c. **Moment vector**

       The moment carried by the bond is accumulated incrementally from the
       relative rotation velocity, with the rolling stiffness :math:`k_r`
       (keyword ``krInnerBond`` or ``krOuterBond``):

       .. math::

          \Delta \underline{M} = k_r \left( \underline{\omega}_j - \underline{\omega}_i \right) \Delta t

       Unlike the resistant moment of a *contact*, this moment is **not**
       capped: as long as the bond exists it grows freely, and it is the
       rupture criterion below that limits it through :math:`M_0`.

    d. **Rupture criterion**

       The rupture criterion determines if the cohesive bond should break. If the yield value is greater than 0, the bond is irreversibly broken (actually, the breakage of the whole interface is scheduled as soon as possible). The yield function reads:

       .. math::

          \varphi = \left( \frac{\Vert \underline{f}_t \Vert}{f_t^0} \right)^p + \left(\frac{\Vert\underline{M}\Vert}{M_0}\right)^p - \frac{f_n}{f_n^0} - 1 

       with the strengths :math:`f_n^0`, :math:`f_t^0`, :math:`M_0` and the
       exponent :math:`p` (keywords ``fn0*Bond``, ``ft0*Bond``, ``mom0*Bond``
       and ``pow*Bond``). The breakage is postponed until all the forces of
       the step are computed, to avoid an asymmetry in the tangential forces.


2. **Contact**

    When the bodies are not bonded (or once their bond is broken), the law is
    that of ``Default``, with the ``*Contact`` parameters and the contact
    partnership, with two differences: the Coulomb threshold uses the total
    normal force, :math:`\Vert \underline{f}_t \Vert \leq \mu f_n`, and the
    moment is computed even when :math:`k_r = 0`.


Bonded Cell Method (keyword ``BCM``)
------------------------------------

This law was introduced to model a cohesive solid discretised into cells that
share *faces*, rather than a packing of grains glued at a few points. The
interface between two cells carries an **area**, and the rupture is governed by
a **fracture energy** rather than by a force threshold.

.. important::

   ``BCM`` only makes sense once the interfaces have been created by the
   pre-processing command ``stickBCM``, which is the command that computes the
   area of each interface (see :ref:`prePro`). All the interfaces built that
   way are flagged as *inner*, so ``BCM`` reads only the ``*InnerBond``
   parameters.

Cohesive bond
^^^^^^^^^^^^^

The values of :math:`k_n` and :math:`k_t` come from ``knInnerBond`` and
``ktInnerBond``, or from the interface itself when ``ParamsInInterfaces`` is 1.
By default, they are the stiffnesses of **each bond**. With
``bondStiffnessPerInterface 1``, they are the stiffnesses of the **whole
interface**, shared between its :math:`N_b` bonds:

.. math::

   k_n \leftarrow \frac{k_n}{N_b} \qquad k_t \leftarrow \frac{k_t}{N_b}

so that refining the discretisation of an interface does not stiffen it. The
damping coefficient of the bonds and the breakage energy below use the same
weighted stiffnesses. For the cubic cells of :ref:`softParticles`, which have
4 bonds per face, ``knInnerBond 0 0 1e5`` per bond and ``knInnerBond 0 0 4e5``
per interface give the same result.

The normal force combines the elastic response, measured from the distance
:math:`d_n^0` of the bond when it was created by ``stickBCM`` (each bond has its
own, so that a bond carries no force in its initial position), with the viscous
damping:

.. math::

   f_n = -k_n \left( d_n - d_n^0 \right) + c\, v_n

and the tangential force is the sum of an elastic part, accumulated
incrementally, and of a viscous part:

.. math::

   \underline{f}_t^e \leftarrow \underline{f}_t^e + k_t\, \underline{v}_t \Delta t
   \qquad
   \underline{f}_t = \underline{f}_t^e + c\, \underline{v}_t

The viscous part of the previous step is removed before the elastic increment
is added, so that the viscosity is not cumulated over time.

Neither the tangential force nor the moment is capped while the bond holds: the
bond is elastic until it breaks.

Rupture criterion
^^^^^^^^^^^^^^^^^

The criterion is energetic and applies to the interface **as a whole**, not to
its individual bonds. At each step, the elastic energy stored in the interface
is summed over its bonds. Only the bonds in *tension* contribute to the normal
part:

.. math::

   E_n = \sum_{b\,\mid\, d_n^b > d_n^{0,b}} \frac{1}{2} k_n \left( d_n^b - d_n^{0,b} \right)^2
   \qquad
   E_t = \sum_{b} \frac{1}{2} k_t \left\Vert \underline{s}_b \right\Vert^2

where :math:`\underline{s}_b` is the accumulated tangential sliding of the bond.
The whole interface breaks as soon as

.. math::

   E_n + E_t > 2\, A\, G_c

with :math:`A` the area of the interface and :math:`G_c` the fracture energy
(keyword ``gcInnerBond``, or ``gcOuterBond`` for an outer interface). The factor
2 accounts for the two surfaces created by the crack.

The stiffnesses and :math:`G_c` are those of the force computation: per bond or
per interface (``bondStiffnessPerInterface``), and taken from the interface
itself when ``ParamsInInterfaces`` is 1.

The area, and the :math:`d_n^0` and :math:`\underline{s}_b` of each bond, are
saved in the conf-files (see :ref:`syntaxConf`), so that a computation restarted
from a conf-file breaks its interfaces at the same time as the uninterrupted
one.

When the criterion is met the interface is not broken immediately: it is added
to a list and all its bonds are released once every force has been computed.
This postponement avoids introducing an asymmetry in the incrementally computed
tangential forces.

Contact
^^^^^^^

Once a bond is broken, or between two cells that were never glued, the
interaction falls back to a law identical to ``Default``: linear elastic normal
repulsion with viscous damping and a positivity clamp, Coulomb friction capped
at :math:`\mu f_n`, and a resistant moment capped at :math:`\mu_r f_n`. The
parameters are then ``knContact``, ``ktContact``, ``muContact``, ``krContact``
and ``murContact``, weighted by the contact partnership when one is active.

The interface damage accumulated over a simulation can be followed with the
``TrackDamage`` data extractor (see :ref:`dataExtractors`).


Slow landslides (keyword ``GeoVisc``)
-------------------------------------

``GeoVisc`` is a variant of ``Default`` designed for slow landslides, where the
relevant resistance is viscous rather than elastic-frictional. Its specific
feature is to make the friction force **viscous** while still yielding on the
Coulomb cone.

Normal component
^^^^^^^^^^^^^^^^

Identical to ``Default``: linear elastic repulsion plus viscous damping, with
the total clamped to remain positive.

.. math::

   f_n = \left[ -k_n d_n + c\, v_n \right]_{\geq 0}

The keywords are ``knContact`` and ``en2Contact``.

Tangential component
^^^^^^^^^^^^^^^^^^^^

This is where the law departs from ``Default``. The tangential force is *not*
accumulated incrementally from a tangential stiffness; it is proportional to the
current sliding velocity, through a tangential viscosity
:math:`\eta_t` (keyword ``viscTContact``):

.. math::

   \underline{f}_t = \eta_t\, \underline{v}_t

It is then projected back onto the Coulomb cone, using the **elastic** part of
the normal force as the reference:

.. math::

   \left\Vert \underline{f}_t \right\Vert \leq \mu \left| f_n^\text{elas} \right|

.. note::

   Using :math:`f_n^\text{elas}` rather than the total :math:`f_n` keeps the
   threshold free of the viscous fluctuations of the normal force. Compiling
   with ``ALLOW_NEGATIVE_THRESHOLDS`` switches the reference to the total
   :math:`f_n`.

Because the tangential force has no memory, ``ktContact`` is not used by this
law: a body at rest carries no tangential force, and the model cannot sustain a
static shear. It is meant for creeping motion, not for equilibrium.
The exercise 4 of the tutorial ``examples/columnCollapse`` shows it: a deposit
that would stop with ``Default`` keeps creeping with ``GeoVisc``.

Moment vector
^^^^^^^^^^^^^

As in ``Default``, and only when :math:`k_r > 0`:

.. math::

   \Delta \underline{M} = \left[ k_r \left( \underline{\omega}_j - \underline{\omega}_i \right) \Delta t \right]_{\pm \mu_r f_n}

The keywords are ``krContact`` and ``murContact``.


Summary of the parameters used by each law
------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 13 13 13 13 14

   * - Parameter
     - ``Default``
     - ``Avalanche``
     - ``StickedLinks``
     - ``BCM``
     - ``GeoVisc``
   * - ``knContact``, ``en2Contact``
     - yes
     - yes
     - yes
     - yes
     - yes
   * - ``ktContact``
     - yes
     - yes
     - yes
     - yes
     - no
   * - ``muContact``, ``krContact``, ``murContact``
     - yes
     - yes
     - yes
     - yes
     - yes
   * - ``viscTContact``
     - no
     - no
     - no
     - no
     - yes
   * - ``kn/ktInnerBond``, ``kn/ktOuterBond``
     - no
     - no
     - yes
     - inner only
     - no
   * - ``krInnerBond``, ``krOuterBond``
     - no
     - no
     - yes
     - no
     - no
   * - ``en2InnerBond``, ``en2OuterBond``
     - no
     - no
     - yes
     - inner only (damping set by ``stickBCM``)
     - no
   * - ``fn0/ft0/mom0/pow*Bond``
     - no
     - no
     - yes
     - no
     - no
   * - ``gcInnerBond``, ``gcOuterBond``
     - no
     - no
     - no
     - yes
     - no

A parameter that a law does not read may be left undefined in the conf-file.
Conversely, every pair of groups that can come into contact must have the
parameters of the selected law defined, otherwise they default to zero.
