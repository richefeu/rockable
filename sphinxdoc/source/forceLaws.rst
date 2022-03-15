.. _Force-laws:

Force laws
==========

We provide here the interaction models existing in ``Rockable``.
The equations are explained with the associated parameters.
These parameters ``[?PARAMETER?]`` are to be defined in the input file in the following way: 

- ``[?PARAMETER?]`` (*int*) **firstGroupNumber** (*int*) **secondGroupNumber** (*double*) **value**


Default (keywork ``Default``)
-----------------------------

This is the force law to be used in most cases. 
In short, the model includes elastic linear contact, normal viscosity, constant normal contact force, 
Coulomb friction, and rolling resistance


Normal component
""""""""""""""""
The elastic part of the normal contact force is

.. math::
   f_n^{el} = -k_n d_n

where :math:`d_n`` is the normal negative or null distance, so :math:`f_n^{el} \geq 0`. The value of :math:`k_n` can be set with the keyword ``knContact``

.. math::
   f_n^{visc} = \alpha_n \sqrt{2 m_{eff}} v_n


Tangential component
""""""""""""""""""""

.. math::
   \Delta f_t^{el} = -k_t (v_t \delta t)


Avalanche (keywork ``Avalanche``)
---------------------------------

TODO


Breakable elastic solid bonds 
-----------------------------

TODO



