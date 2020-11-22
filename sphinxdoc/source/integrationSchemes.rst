Integration schemes
===================

In a ``conf-file``, the integration scheme to be used is set with the keyword ``Integrator`` followed by the name of the scheme (defined below). The position, velocity and acceleration vectors is noted :math:`\vec{r}`, :math:`\vec{v}` and :math:`\vec{a}`, respectively; and for rotations they are :math:`\mathbf{Q}`, :math:`\vec{\omega}` and :math:`\vec{\dot{\omega}}`. 

For all these schemes, the driven bodies with force or moment are updated similarly to the chosen scheme, and the driven bodies with velocity imposed (translation or rotation) are updated that way:

.. math::
   \begin{cases}
   \vec{r}(t+\delta t) &\leftarrow \vec{r}(t) + \vec{v}^\mbox{imp}(t) \delta t \\
   \mathbf{Q}(t+\delta t) &\leftarrow \mathbf{Q}(t) + \mathbf{\dot{Q}}\left(\vec{\omega}^\mbox{imp}(t)\right) \delta t 
   \end{cases}

Euler scheme (``Euler``)
------------------------

This is the most simple scheme:

.. math::
   \begin{cases}
   \vec{r}(t+\delta t) &\leftarrow \vec{r}(t) + \vec{v}(t) \delta t \\
   \vec{v}(t+\delta t) &\leftarrow \vec{v}(t) + \vec{a}(t) \delta t \\
   \vec{a}(t+\delta t) &\leftarrow \ldots
   \end{cases}
   
And the same for rotations

.. math::
   \begin{cases}
   \mathbf{Q}(t+\delta t) &\leftarrow \mathbf{Q}(t) + \mathbf{\dot{Q}}\left(\vec{\omega}(t)\right) \delta t \\
   \vec{\omega}(t+\delta t) &\leftarrow \vec{\omega}(t) + \vec{\dot{\omega}}(t) \delta t \\
   \vec{\dot{\omega}}(t+\delta t) &\leftarrow \ldots
   \end{cases}

Velocity-Verlet scheme (``velocityVerlet``, Default)
----------------------------------------------------

.. math::
   \begin{cases}
   \vec{r}(t+\delta t) &\leftarrow \vec{r}(t) + \vec{v}(t) \delta t + \vec{a}(t) \frac{\delta t^2}{2} \\
   \vec{v}(t+\frac{\delta t}{2}) &\leftarrow \vec{v}(t) + \vec{a}(t) \delta t \\
   \vec{a}(t+\delta t) &\leftarrow \ldots \\
   \vec{v}(t+\delta t) &\leftarrow \vec{v}(t+\frac{\delta t}{2}) + \vec{a}(t+\frac{\delta t}{2}) \delta t \\
   \end{cases}

.. math::
   \begin{cases}
   \mathbf{Q}(t+\delta t) &\leftarrow \mathbf{Q}(t) + \mathbf{\dot{Q}}\left(\vec{\omega}(t)\right) \delta t + \mathbf{\ddot{Q}}\left(\vec{\omega}(t),\vec{\dot{\omega}}(t)\right) \frac{\delta t^2}{2} \\
   \vec{\omega}(t+\frac{\delta t}{2}) &\leftarrow \vec{\omega}(t) + \vec{\dot{\omega}}(t) \delta t \\
   \vec{\dot{\omega}}(t+\delta t) &\leftarrow \ldots \\
   \vec{\omega}(t+\delta t) &\leftarrow \vec{\omega}(t+\frac{\delta t}{2}) + \vec{\dot{\omega}}(t+\frac{\delta t}{2}) \delta t \\
   \end{cases}


Beeman scheme (``Beeman``)
--------------------------



Runge-Kutta-Nyström 4th-order scheme (``RungeKutta4``)
------------------------------------------------------


