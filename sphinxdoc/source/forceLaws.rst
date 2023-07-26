.. _Force-laws:

Force laws
==========

We provide here the interaction models existing in ``Rockable``.
The equations are explained with the associated parameters.
The parameters ``[PARAMETER]`` are to be defined in the input conf-file in the following way: 

- ``[PARAMETER]`` (*int*) **firstGroupNumber** (*int*) **secondGroupNumber** (*double*) **value**


Default model (keywork ``Default``)
-----------------------------------

This is the force law to be used in most cases. 
In short, the model includes elastic linear normal force for contact, normal viscosity force, constant normal cohesion force at contact, 
Coulomb friction tangent force, and rolling resistance moment


Normal component
""""""""""""""""
The elastic part of the normal contact force is

.. math::

   f_n^{el} = -k_n d_n

where :math:`d_n \leq 0` is the normal distance, so :math:`f_n^{el} \geq 0`.
When :math:`dn > 0`, this force is zero.
The value of :math:`k_n` can be set with the keyword ``knContact``

The viscuous part of the normal contact force is

.. math::

   f_n^{visc} = \alpha_n \sqrt{2 m_\text{eff}} v_n

where :math:`v_n` is the relative velocity,  :math:`m_\text{eff}=(m_i m_j)/(m_i+m_j)` is the effective mass, 
and :math:`\alpha_n \in [0, 1[` is the rate of normal viscuous damping. There are two ways to set the value of :math:`\alpha_n`:

1. The keyword ``en2Contact`` that is the energy normal restitution rate (:math:`e_n^2`). In this case, the viscuous damping rate will be set associated

.. math::

   \alpha_n = \frac{- \ln e_n}{\sqrt{\ln^2 e_n + \pi^2}}

2. The keyword ``en2ContactFromViscRate`` that will  

.. math:: 

   e_n^2 = \exp \left(-\frac{\alpha_n \pi}{\sqrt{1 - \alpha_n^2}}\right)

Because of the viscuous part of the normal force, the total normal force :math:`f_n = f_n^{el} + f_n^{visc}` can be negative. 
Although the physical meaning of this negative value is arguable, 
the normal force is by default restricted to remain positive or zero.

.. note:: After discussion with Farhang, an option will be added to let force not to be restricted 
          (let it be negative if it wants)


Tangential component
""""""""""""""""""""

.. TODO

.. math::
   \Delta f_t^{el} = -k_t (v_t \delta t)


Moment  component
"""""""""""""""""

.. TODO


Law for rock avalanches (keywork ``Avalanche``)
-----------------------------------------------

This is historically the first law that has been implemented in ``Rockable`` (actually in ``DEMbox``, its ancestor).

.. TODO


Breakable elastic solid bonds (keywork ``StickedLinks``)
--------------------------------------------------------

This law has been introduced since the PhD work of Marta Stasiak.

.. TODO



