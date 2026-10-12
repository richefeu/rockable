.. _Dissipation:

Global dissipation features
===========================

Cundall damping
---------------

In DEM, Cundall damping is a numerical technique employed to introduce artificial dissipation and stabilize quasistatic simulations of granular materials. Unlike physical processes that naturally dissipate energy, Cundall damping is a purely numerical approach aimed at controlling the excessive accumulation of energy and preventing unrealistic behavior in the simulation.

Cundall damping does not add a force: it rescales the resultant force (and
moment) of each body, reducing it when it drives the motion and increasing it
when it opposes the motion. The bodies are therefore slowed down without any
reference to a physical dissipative process, which helps reach a quasi-static
equilibrium.

In ``Rockable``, both translations and rotations are affected, for the free
bodies and for the force-driven components of the driven bodies. By default, the
test is made on the whole vectors (:math:`\underline{F} \cdot \underline{v}`);
if the code has been compiled with the ``COMPONENTWISE_NUM_DAMPING`` option, each
component is processed separately. It is not clear at all what the right
solution should be. To activate this damping, the keyword is
``numericalDampingCoeff``, set with a value between 0 and 1. Here is how the
damping works:

.. code-block:: c++
   :caption: Pseudo-code
	 
   if (force * velocity > 0.0) {
	   force *= (1.0 - numericalDampingCoeff)
   } else {
	   force *= (1.0 + numericalDampingCoeff)
   }

.. note::

   With a periodic cell (see :ref:`periodicBoundaryConditions`), the sign is
   tested with the velocity of each body relative to the mean deformation of the
   cell, so that the damping does not act on the affine motion imposed by the
   cell.

There is no universally fixed or commonly used value for the Cundall coefficient since it depends on factors such as the material properties, particle interactions, and the desired behavior of the simulation. The appropriate value of ``numericalDampingCoeff`` is typically determined through calibration and validation against experimental data or known behavior of the material being simulated. It is common practice to start with a small value and gradually increase it until the desired level of damping is achieved. The specific value used can vary widely, ranging from very small values (*e.g.*, 0.001) to larger values (*e.g.*, 0.1) or even higher, depending on the simulation requirements and the desired damping effect.

It's important to note that the choice of the Cundall damping coefficient involves a trade-off. A larger value can provide stronger damping but may also introduce more artificial dissipation and affect the overall accuracy of the simulation. Therefore, it is often necessary to conduct sensitivity analyses and validation studies to determine the appropriate value of this coefficient for a specific DEM simulation.
 

Velocity barrier
----------------

An alternative solution to mitigate the excessive free motion of particles in DEM simulations is to apply a barrier function to translational and rotational velocities. This approach aims to impose limits on the particle velocities, preventing them from exceeding certain thresholds and thus restraining overly fast or unstable movements. This is an original idea of Farhang Radjai. 

The barrier function is typically defined as a threshold function that specifies a maximum allowable value for the translational and rotational velocities. If a particle surpasses these limits, its velocity is appropriately reduced to adhere to the imposed constraints. Here is the pseudo-code for each component of force (or equivalently, of acceleration):

.. code-block:: c++
   :caption: Pseudo-code

    ratio = pow(fabs(velocity / Barrier), Exponent);
    force *= (1.0 - ratio) / (1.0 + ratio);


The advantage of this method is that it can be used to specifically control particle velocities without introducing artificial dissipation forces. Hence, it provides an alternative to dissipate energy, in particular for deposition of particles. The keywords are: ``VelocityBarrier`` with ``VelocityBarrierExponent`` for translation, and ``AngularVelocityBarrier`` with ``AngularVelocityBarrierExponent`` for rotations. A barrier is active when its value is strictly positive.

In the source, the factor is applied to the **accelerations**, component by
component. It is close to 1 for a velocity much smaller than the barrier, zero
at the barrier, and negative beyond it (it tends to -1), so that a velocity
exceeding the barrier is brought back towards it.


Particles in viscous fluid
--------------------------

The presence of a surrounding fluid can be represented by a drag force applied
to each body, which is the ``BodyForce`` named ``ViscousFluid``, set with the
density of the fluid:

.. code-block:: text
   :caption: input.txt

   BodyForce ViscousFluid 1000

The intended force is a quadratic drag, :math:`\frac{1}{2} C_X \rho V^{2/3} v^2`
with :math:`C_X = 0.47` as for a sphere. In its current implementation, it does
not oppose the motion: see the warning in :ref:`bodyForces` before using it.
