Integration schemes
===================

In a ``conf-file``, the integration scheme to be used is set with the keyword ``Integrator`` followed by the name of scheme (defined above).

``Euler``
---------

This is the most simple scheme:

.. math::
   \begin{cases}
   \vec{r}(t+\delta t) &= \vec{r}(t) + \vec{v}(t) \delta t \\
   \vec{v}(t+\delta t) &= \vec{v}(t) + \vec{a}(t) \delta t \\
   \vec{a}(t+\delta t) &= \mbox{resultantForce}(\ldots)
   \end{cases}

``velocityVerlet`` (Default)
----------------------------

.. math::
   \begin{cases}
   \vec{r}(t+\delta t) &= \vec{r}(t) + \vec{v}(t) \delta t + \vec{a}(t) \frac{\delta t^2}{2} \\
   \vec{v}(t+\delta t) &= \vec{v}(t) + \vec{a}(t) \delta t
   \end{cases}


