.. _Dissipation:

Global dissipation features
===========================

Cundall damping
---------------

In DEM, Cundall damping is a numerical technique employed to introduce artificial dissipation and stabilize quasistatic simulations of granular materials. Unlike physical processes that naturally dissipate energy, Cundall damping is a purely numerical approach aimed at controlling the excessive accumulation of energy and preventing unrealistic behavior in the simulation.

By incorporating Cundall damping, a damping force proportional to the relative velocity between particles is applied, which leads to the gradual reduction of particle velocities over time. This numerical dissipation helps stabilize the simulation, particularly in situations where quasi-static equilibrium is desired.

To summarize, Cundall damping in DEM is a numerical dissipation technique implemented to stabilize quasi-static simulations of granular materials, by artificially introducing a damping force that reduces particle velocities without directly modeling physical dissipative processes.

In ``Rockable``, both translations and rotations are affected. If the code has been compiled with ``COMPONENTWISE_NUM_DAMPING`` option, the damping of each component will be processed separately. It is not clear at all what the right solution should be. To activate this damping, the keyword is ``numericalDampingCoeff``, set with a positive value that can vary depending on the specific simulation and the desired level of damping. Basically, here is what the damping works:

.. code-block:: c++
   :caption: Pseudo-code
	 
   if (force * velocity > 0.0) {
	   force *= (1.0 - numericalDampingCoeff)
   } else {
	   force *= (1.0 + numericalDampingCoeff)
   }

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


The advantage of this method is that it can be used to specifically control particle velocities without introducing artificial dissipation forces. Hence, it provides an alternative to dissipate energy, in particular for deposition of particles. The keywords are: ``velocityBarrier`` with ``velocityBarrierExponent`` for translation, and ``angularVelocityBarrier`` with ``angularVelocityBarrierExponent`` for rotations.


Particles in viscous fluid
--------------------------

This dissipation strategy involves adding a viscosity term to the particle movements in order to simulate the presence of a viscous fluid surrounding the particles. This approach is commonly used to model interactions between particles and a surrounding fluid, such as in fluid-structure interaction simulations. In this strategy, a damping force is introduced that is proportional to the velocity of the particles, similar to the behavior of a viscous fluid. The damping force opposes the particle's motion and causes energy dissipation over time.

Mathematically, the damping force can be expressed as:

.. code-block:: c++
   :caption: Pseudo-code
	 
   F_damping = -viscuous_damping * velocity


By incorporating this damping force, the particle's motion is gradually slowed down, and its kinetic energy is dissipated. This mimics the effect of a surrounding viscous fluid, which would resist and dampen the particle's motion over time.

This is actually implemented as the ``BodyForce`` named ``ViscousFluid`` in ``Rockable``.

To set a more or less appropriate value of the viscuous damping coefficient, the user can set a fluid density :math:`\rho` that way:

.. code-block:: text
   :caption: input.txt
	 
   ViscousFluid 1000

And the viscuous damping value will be set to :math:`\frac{1}{2} C_X\ \rho\ V^{2/3}`, with :math:`C_X = 0.47` just like if the particle was a sphere.




