.. _bodyForces:

Body forces
===========

A *body force* acts on each body individually, from its own state, without
involving any interaction. Gravity is the obvious example, but it is handled
separately by the ``gravity`` keyword; the ``BodyForce`` mechanism covers the
other cases.

.. code-block:: text
   :caption: input.txt

   BodyForce ViscousFluid 1000

.. important::

   There can be only **one** body force at a time. Declaring a second
   ``BodyForce`` replaces the first. Unlike the data extractors, the active body
   force *is* written back into the conf-files, so it survives a restart.


``ViscousFluid`` (*double*) **fluidDensity**
--------------------------------------------

Drag exerted by a surrounding fluid of density :math:`\rho`. The intended model
is the usual quadratic drag, applied component by component:

.. math::

   \underline{F} = \frac{1}{2} C_X\, \rho\, V^{2/3}\, \underline{v}^{\,2}

where :math:`V` is the volume of the body (homothety included) and
:math:`C_X = 0.47` is the drag coefficient of a sphere. The exponent
:math:`2/3` turns the volume into a frontal area, so that a body twice as large
in every direction feels four times the drag.

This is the way to deposit a sample in a liquid, or simply to dissipate energy
in a physically motivated manner rather than with the numerical damping
described in :ref:`Dissipation`.

.. warning::

   In the current implementation the force is built component by component from
   the **squared** velocity, without restoring the sign of the velocity:

   .. code-block:: c++

      force = 0.5 * C_X * rho * pSurface * vec3r(vx*vx, vy*vy, vz*vz);

   Each component is therefore always positive, and the force does not oppose
   the motion as a drag should: a body moving along :math:`-x` is pushed
   further along :math:`+x`. The source itself carries a note saying that this
   solution has never been tested. Treat this body force as experimental, and
   check its effect on your case before using it.


``AttractingPoint`` (*vec3r*) **point** (*double*) **acceleration**
-------------------------------------------------------------------

A constant acceleration pulling every body towards a fixed point of space:

.. math::

   \underline{F} = m\, a\, \frac{\underline{x}_\text{point} - \underline{x}}{\Vert \underline{x}_\text{point} - \underline{x} \Vert}

The magnitude does not decrease with distance: this is a *centripetal* field of
constant intensity, not a gravitational one. It is convenient to compact a
sample towards a centre without any wall, or to hold a heap on a curved surface.

The three first values are the coordinates of the point, the fourth is the
acceleration.


``PreferredDirection`` (*vec3r*) **axisBody** (*vec3r*) **axis** (*double*) **momentMax**
------------------------------------------------------------------------------------------

A restoring moment that brings a given axis of the body onto a fixed direction
of space. It is used to model bodies that tend to align, such as elongated
grains in a flow, or to keep a body upright.

``axisBody`` is the axis carried by the body, expressed in its **own frame**;
``axis`` is the target direction, in the global frame. Both are normalised
automatically. The moment is proportional to the misalignment angle
:math:`\theta`:

.. math::

   \underline{M} = k_r\, \theta\, \underline{u}
   \qquad\text{with}\qquad
   k_r = \frac{2 M_\text{max}}{\pi}

where :math:`\underline{u}` is the unit vector of the rotation bringing one axis
onto the other. The stiffness is set so that the moment reaches
:math:`M_\text{max}` when the body axis is perpendicular
(:math:`\theta = \pi/2`) to the target direction, which makes ``momentMax`` a
directly interpretable parameter.

.. note::

   The target direction is treated as a **director, not as a vector**: when the
   body axis points away from ``axis``, the angle is folded back through
   :math:`\pi - \theta` and the sign of the moment is flipped. A body therefore
   aligns with whichever of :math:`+\underline{a}` or :math:`-\underline{a}` is
   the closer, and :math:`\theta` never exceeds :math:`\pi/2`. Reversing the
   sign of ``axis`` in the conf-file changes nothing.


Relation with the other dissipation mechanisms
----------------------------------------------

``ViscousFluid`` is one of three ways of removing energy from a sample, and the
choice matters:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Mechanism
     - Nature
   * - ``numericalDampingCoeff``
     - Purely numerical (Cundall). Acts on the forces, at every contact and on
       every body. Fast to bring a sample to equilibrium, but has no physical
       meaning.
   * - ``VelocityBarrier``
     - Caps the velocities without adding any dissipative force. Useful for a
       deposition, leaves the contact law untouched.
   * - ``BodyForce ViscousFluid``
     - A real drag force, with a physical parameter. The only one of the three
       that can be defended in a published result as part of the model.

They are described together in :ref:`Dissipation`.
