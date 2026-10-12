.. _IntegrationSchemes:

Integration schemes
===================

In a ``conf-file``, the integration scheme to be used is set with the keyword ``Integrator`` 
followed by the name of the scheme (defined below). The position, velocity and acceleration 
vectors is noted :math:`\underline{x}`, :math:`\underline{v}` and :math:`\underline{a}`, respectively; 
and for rotations they are :math:`\mathbf{Q}` (a quaternion), :math:`\underline{\omega}` and :math:`\underline{\dot{\omega}}`. 

For all these schemes, the driven bodies with force or moment are updated similarly to the chosen scheme, 
and the driven bodies with velocity imposed (translation or rotation) are updated that way:

.. math::

   \begin{cases}
   \underline{x}_{t+\Delta t} &= \underline{x}_{t} + \underline{v}^\mbox{imposed}_{t} \Delta t \\
   \mathbf{Q}_{t+\Delta t} &= \mathbf{Q}_{t} + \mathbf{\dot{Q}}\left(\underline{\omega}^\mbox{imposed}_{t}\right) \Delta t 
   \end{cases}

Euler scheme (keyword ``Euler``)
--------------------------------

The Euler scheme is a simple and intuitive time integration method, but it has limitations. One of the main limitations is that it is conditionally stable. Therefore, the time step size must be carefully chosen to ensure numerical stability. Additionally, the explicit Euler method is only first-order accurate, meaning that errors in position and velocity tend to accumulate over time, especially for stiff systems. This can lead to inaccuracies in the simulation, especially for long integration times or complex systems.


For a discrete particle in DEM, the equations used to update its position and velocity using the Euler scheme are as follows:


**Position update**


The new position and orientation of a particle, :math:`\underline{x}_{t+\Delta t}` and :math:`\mathbf{Q}_{t+\Delta t}`, at the next time step (:math:`t+\Delta t`) are obtained from its current position and orientation, :math:`\underline{x}_t` and :math:`\mathbf{Q}_{t}`, and its velocity, :math:`\underline{v}_t` and :math:`\mathbf{\dot{Q}}_t`, at the current time step (:math:`t`). The position update equation is given by:

.. math::

   \begin{cases}
   \underline{x}_{t+\Delta t} &= \underline{x}_t + \underline{v}_t \Delta t \\
	 \mathbf{Q}_{t+\Delta t}    &= \mathbf{Q}_t + \mathbf{\dot{Q}}(\underline{\omega}_{t}) \Delta t
	 \end{cases}
	 
where :math:`\Delta t` is the time step size, which is a user-defined parameter.


**Velocity update**


The new velocity of a particle, :math:`\underline{v}_{t+\Delta t}` and :math:`\underline{\omega}_{t+\Delta t}`, at the next time step (:math:`t+\Delta t`) is obtained from its current velocity, :math:`\underline{v}_t` and :math:`\underline{\omega}_t`, and the acceleration, :math:`\underline{a}_t` and :math:`\underline{\dot{\omega}}_t`, at the current time step (:math:`t`). The velocity update equation is given by:

.. math::

   \begin{cases}
   \underline{v}_{t+\Delta t}      &= \underline{v}_t + \underline{a}_t \Delta t \\
	 \underline{\omega}_{t+\Delta t} &= \underline{\omega}_t + \underline{\dot{\omega}}_t \Delta t
	 \end{cases}
	 
where :math:`\underline{a}_t` (and :math:`\underline{\dot{\omega}}_t`) is the acceleration calculated from the forces (and moments) acting on the particle at time step :math:`t`.



Velocity-Verlet scheme (keyword ``velocityVerlet``, Default)
------------------------------------------------------------


The Velocity Verlet method has several advantages. It is second-order accurate, meaning that it reduces the numerical errors in position and velocity and improves the overall accuracy of the simulation. Additionally, the Velocity Verlet method is time-reversible, which means that the integration process can be reversed without introducing significant errors. Moreover, it conserves energy, making it particularly suitable for conservative systems like those encountered in DEM.

.. warning:: Because this scheme is the default one, most new features may be implemented here and maybe not in other                schemes. Check the code to be sure.

For a discrete particle in DEM, the equations used to update its position and velocity using the Velocity Verlet method are as follows:


**Position update**


The new position of a particle, :math:`\underline{x}_{t+\Delta t}` and :math:`\mathbf{Q}_{t+\Delta t}`, at the next time step (:math:`t+\Delta t`) is obtained from its current position, :math:`\underline{x}_t` and :math:`\mathbf{Q}_{t}`, its velocity, :math:`\underline{v}_t` and :math:`\underline{\omega}_t`, and its acceleration, :math:`\underline{a}_t` and :math:`\underline{\dot{\omega}}_t`, at the current time step (:math:`t`). The position update equation is given by:

.. math::

   \begin{cases}
   \underline{x}_{t+\Delta t} &= \underline{x}_t + \underline{v}_t \Delta t + \frac{1}{2}\underline{a}_t \Delta t^2\\
	 \mathbf{Q}_{t+\Delta t}    &= \mathbf{Q}(t) + \mathbf{\dot{Q}}\left(\underline{\omega}_t\right) \Delta t + \frac{1}{2} \mathbf{\ddot{Q}}\left(\underline{\omega}_t,\underline{\dot{\omega}}_t\right) \Delta t^2 \\
   \end{cases}
	 
where :math:`\Delta t` is the time step size, which is a user-defined parameter.


**Velocity update**

The new velocity of a particle, :math:`\underline{v}_{t+\Delta t}` and :math:`\underline{\omega}_{t+\Delta t}`, at the next time step (:math:`t+\Delta t`) is obtained from its current velocity, :math:`\underline{v}_t` and :math:`\underline{\omega}_t`, and an intermediate acceleration between the current and the next time steps. The velocity update equation is given by:

.. math::

   \begin{cases}
   \underline{v}_{t+\Delta t} &= \underline{v}_t + \frac{1}{2} \left( \underline{a}_t + \underline{a}_{t+\Delta t} \right) \Delta t \\
   \underline{\omega}_{t+\Delta t} &= \underline{\omega}(t) + \frac{1}{2} \left( \underline{\dot{\omega}}_t + \underline{\dot{\omega}}_{t+\Delta t} \right) \Delta t \\
	 \end{cases}

where :math:`\underline{a}_{t+\Delta t}` (and :math:`\underline{\dot{\omega}}_{t+\Delta t}`)  is the acceleration calculated from the forces (and moments) acting on the particle at time step :math:`t+\Delta t`.


Beeman scheme (keyword ``Beeman``)
----------------------------------


The Beeman method is a second-order accurate integration technique, which means it provides better accuracy compared to the first-order explicit Euler (or Velocity Verlet) method. This leads to more accurate and reliable simulations of particle dynamics in DEM. Additionally, the Beeman method is relatively easy to implement and computationally efficient compared to some higher-order methods.

For a discrete particle in DEM, the equations used to update its position and velocity using the Beeman method are as follows:


In Rockable, it is used as a predictor-corrector, written here for the
translations (the rotations follow the same logic with quaternions).

**Position update**

The new position uses the current and the previous accelerations:

.. math::

   \underline{x}_{t+\Delta t} = \underline{x}_t + \Delta t \, \underline{v}_t + \frac{\Delta t^2}{6} \left( 4 \underline{a}_t - \underline{a}_{t-\Delta t} \right)

**Predicted velocity**

The forces at :math:`t+\Delta t` are computed with a predicted velocity (it
matters for the viscous forces):

.. math::

   \underline{v}^\text{pred}_{t+\Delta t} = \underline{v}_t + \Delta t \left( \frac{3}{2} \underline{a}_t - \frac{1}{2} \underline{a}_{t-\Delta t} \right)

**Corrected velocity**

Once the acceleration :math:`\underline{a}_{t+\Delta t}` is known, the velocity
is corrected:

.. math::

   \underline{v}_{t+\Delta t} = \underline{v}_t + \frac{\Delta t}{12} \left( 5 \underline{a}_{t+\Delta t} + 8 \underline{a}_t - \underline{a}_{t-\Delta t} \right)

The scheme computes the accelerations twice per time step.


Runge-Kutta-Nyström 4th-order scheme (keyword ``RungeKutta4``)
--------------------------------------------------------------

	
The Runge-Kutta-Nyström (RKN) 4th-order scheme is a higher-order accurate integration technique, providing better accuracy compared to lower-order methods like the explicit Euler or Velocity Verlet. This can lead to more accurate and reliable simulations of particle dynamics in DEM. Additionally, the RKN 4th-order scheme has good stability properties and can handle stiff systems more effectively than lower-order methods.


.. warning:: For this scheme, the equations for rotations are not shown but follow the same logic as for other schemes.

The accelerations are computed four times per time step. With
:math:`\underline{k}_n` the accelerations computed at the stage :math:`n`:

.. math::

   \begin{aligned}
   \underline{k}_1 &= \underline{a}\left( \underline{x}_t, \; \underline{v}_t \right) \\
   \underline{k}_2 &= \underline{a}\left( \underline{x}_t + \frac{\Delta t}{2} \underline{v}_t + \frac{\Delta t^2}{8} \underline{k}_1, \; \underline{v}_t + \frac{\Delta t}{2} \underline{k}_1 \right) \\
   \underline{k}_3 &= \underline{a}\left( \underline{x}_t + \frac{\Delta t}{2} \underline{v}_t + \frac{\Delta t^2}{8} \underline{k}_1, \; \underline{v}_t + \frac{\Delta t}{2} \underline{k}_2 \right) \\
   \underline{k}_4 &= \underline{a}\left( \underline{x}_t + \Delta t \, \underline{v}_t + \frac{\Delta t^2}{2} \underline{k}_3, \; \underline{v}_t + \Delta t \, \underline{k}_3 \right)
   \end{aligned}

and the position and the velocity at the end of the step are

.. math::

   \begin{aligned}
   \underline{x}_{t+\Delta t} &= \underline{x}_t + \Delta t \, \underline{v}_t + \frac{\Delta t^2}{6} \left( \underline{k}_1 + \underline{k}_2 + \underline{k}_3 \right) \\
   \underline{v}_{t+\Delta t} &= \underline{v}_t + \frac{\Delta t}{6} \left( \underline{k}_1 + 2 \underline{k}_2 + 2 \underline{k}_3 + \underline{k}_4 \right)
   \end{aligned}

.. note::

   The schemes ``Beeman`` and ``RungeKutta4`` compute the forces more than once
   per time step. The force laws that accumulate a quantity at each
   computation (the tangential forces of the contacts, for instance) then
   accumulate it as many times, and the viscoelastic soft particles (see
   :ref:`softParticles`) are only correct with ``Euler`` and ``velocityVerlet``.
