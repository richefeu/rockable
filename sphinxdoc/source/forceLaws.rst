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


Default model (keywork ``Default``)
-----------------------------------

This is the force law to be used in most cases. 
In short, the model includes elastic linear normal force for contact, normal viscosity force, constant normal cohesion force at contact, 
Coulomb friction tangent force, and rolling resistance moment


1. **Normal component**

  a. The elastic part of the normal contact force is

  .. math::

     f_n^\text{elas} = -k_n d_n

  where :math:`d_n \leq 0` is the normal distance, so :math:`f_n^\text{elas} \geq 0`. When :math:`d_n > 0`, this force is zero. The value of :math:`k_n` can be set with the keyword ``knContact``

  b. The viscuous part of the normal contact force is

  .. math::

     f_n^\text{visc} = \alpha_n \sqrt{2 m_\text{eff}} v_n

  where :math:`v_n` is the relative normal velocity, :math:`m_\text{eff}=(m_i m_j)/(m_i+m_j)` is the effective mass, and :math:`\alpha_n \in [0, 1[` is the rate of normal viscuous damping. There are two ways to set the value of :math:`\alpha_n`:

    i. The keyword ``en2Contact`` that is the energy normal restitution rate (:math:`e_n^2`). In this case, the viscuous damping rate will be set associated

    .. math::

       \alpha_n = \frac{- \ln e_n}{\sqrt{\ln^2 e_n + \pi^2}}

    ii. The keyword ``en2ContactFromViscRate`` that will  

    .. math:: 

       e_n^2 = \exp \left(-\frac{\alpha_n \pi}{\sqrt{1 - \alpha_n^2}}\right)


  .. note:: 

     Because of the viscuous part of the normal force, the total normal force :math:`f_n = f_n^\text{elas} + f_n^\text{visc}` can be negative. Although the physical meaning of this negative value is arguable, the normal force is by default restricted to remain positive or zero.


2. **Tangential component**

The tangential force represents the frictional resistance between particles when they slide against each other. It is calculated based on the relative tangential velocity (:math:`v_t`) and a tangential stiffness parameter (:math:`k_t`, keyword ``ktContact``). The Coulomb friction model is used to limit the tangential force, ensuring that it does not exceed the product of the friction coefficient (:math:`\mu`, keyword ``muContact``) and the normal elastic force (:math:`f_n^\text{elas}`). The tangential force is incrementally updated with the following relation, at each time step :math:`\Delta t`, since the onset of contact:

.. math::

   \Delta f_t = \left [ k_t v_t \Delta t \right ]_{\pm \mu f_n^\text{elas}}

The friction force :math:`f_t` is cancelled as soon as contact is lost. 


3. **Moment  vector**

The resistant moment (:math:`\underline{M}`) accounts for rotational resistance between particles. It is calculated based on the rotational velocity difference (:math:`\underline{\omega}_i - \underline{\omega}_j`) and the rotational stiffness parameter (:math:`k_r`, keyword ``krContact``). Additionally, the moment is limited by the product of the moment length coefficient (:math:`\mu_r`, keyword ``murContact``), the normal force (:math:`f_n`), and the length of the branch vector (:math:`\ell`). The resistant moment is updated, similarly to the friction force, as follows:

.. math::

   \Delta \underline{M} = \left [ k_r (\underline{\omega}_i - \underline{\omega}_j) \Delta t \right ]_{\pm \mu_r \ell f_n^\text{elas}}
   

.. note:: The moment resistance is applied only when :math:`k_r > 0`.
	 



Law for rock avalanches (keywork ``Avalanche``)
-----------------------------------------------

This is historically the first law that has been implemented in ``Rockable`` (actually in ``DEMbox``, its ancestor).

The force-law implemented is specifically designed for simulating rock avalanches at Laboratoire 3SR. It governs the interactions between particles or bodies in a granular material, considering both normal and tangential forces as well as resistant moments due to rotational motion. 

1. **Normal component**

   The normal force (`I.fn`) between particles is calculated to represent the contact forces that resist compression and prevent interpenetration. The calculation involves considering the change in the normal displacement (`I.dn`) and the normal stiffness parameter (`kn`). The normal force may also undergo unloading or loading depending on the change in normal displacement.

2. **Tangential component**

   The tangential force (`I.ft`) represents the frictional resistance between particles when they slide against each other. The calculation considers the tangential velocity (`vt`) and a tangential stiffness parameter (`kt`). Frictional forces are limited by the friction coefficient (`mu`) and the normal force (`I.fn`).

3. **Resistant moment**

   The resistant moment (`I.mom`) accounts for the rotational resistance between particles. The calculation considers the rotational velocity difference (`box->Particles[I.j].vrot - box->Particles[I.i].vrot`) and the rotational stiffness parameter (`kr`). Additionally, there is a correction applied to the branch vector (`branch`) to account for different scenarios when the free body is either particle `i` or particle `j`.

.. note:: Weighted Interaction Parameters. The force-law allows for adjusting the interaction parameters (`kn`, `kt`, `kr`) based on a weight factor (`w`). This weighting is determined by the function `box->ctcPartnership.getWeight`, which can influence the forces and moments between the particles.

.. note:: The code includes conditional blocks for certain options (`FT_CORR`) which may affect the actual computations. Additionally, the equations use specific parameters (`en2` and `mur`), which should be defined and obtained from relevant data sources in the simulation setup.



Breakable elastic solid bonds (keywork ``StickedLinks``)
--------------------------------------------------------

This law has been introduced since the PhD work of *Marta Stasiak*.

Two states can exist: cohesive bond and contact.

1. **Cohesive bond**

    When a cohesive bond exists between two particles, the force-law employs the following equations:

    a. **Normal component**

       The normal force (:math:`f_n = f_n^\text{elas} + f_n^\text{visc}`) between the particles is calculated using
       Hooke's elastic contact, combined with viscous damping:

       .. math::

          f_n^\text{elas} = -k_n (d_n - d_n^0)
        
       where :math:`d_n^0` is the initial normal distance. This time, the parameter :math:`k_n` is set with the keyword ``knInnerBond`` or ``knOuterBond``.
       The viscous normal force is:

       .. math::

          f_n^\text{visc} = \alpha_n \sqrt{2 m_\text{eff}} v_n	

       where :math:`\alpha_n` is indirectly set with the keywork ``en2InnerBond`` or ``en2OuterBond``.

    b. **Tangential component**

       The tangential elastic force (:math:`f_t`) is incrementally calculated based on the relative tangential velocity (:math:`v_t`) and a tangential stiffness parameter (:math:`k_t`, keyword ``ktInnerBond`` or ``ktOuterBond``):

       .. math::

           \Delta f_t = k_t (v_t \Delta t)

    Optionally, an additional term for tangential viscosity (:math:`f_t^\text{visc} = \alpha_n \sqrt{2 m_\text{eff}} v_t`) can be added to the elastic part of the tangential force. However, this approach accumulates viscosity over time and may not be accurate:


    c. **Moment vector**
    
    TODO

    d. **Rupture criterion**

       The rupture criterion determines if the cohesive bond should break. If the yield value is greater than 0, the bond is irreversibly broken (actually, the breakage of the whole interface is scheduled as soon as possible). The yield function reads:

       .. math::

          \varphi = \left( \frac{\vert f_t \vert}{f_t^0} \right)^p + \left(\frac{\Vert\underline{M}\Vert}{M_0}\right)^p - \frac{f_n}{f_n^0} - 1 


2. **Contact**

    When particles are in contact (no cohesive bond), the force-law employs the following equations:

    a. **Normal component**

       The normal force (`I.fn`) between the particles is calculated similar to the cohesive bond case, but with different interaction parameters (`kn`, `kt`, `damp`, etc.):

       .. math::

           \text{vn} = \mathbf{I.vel} \cdot \mathbf{I.n}

           \text{fne} = -\text{kn} \cdot \text{I.dn}

           \text{fnv} = \text{damp} \cdot \text{vn}

           \text{I.fn} = \text{fne} + \text{fnv}

    b. **Tangential component** (Friction)

       The tangential force (`I.ft`) due to friction between the particles is calculated in a similar manner to the cohesive bond case but with different interaction parameters (`kt`, `mu`, etc.).

    c. **Resistant moment vector**

       If the contact allows for rotational motion (`kr` parameter), a resistant moment (`I.mom`) is calculated based on the relative rotational velocity of the particles.

Note: The code includes conditional blocks for certain options (`FT_CORR` and `box->ctcPartnership.getWeight != nullptr`) which may affect the actual computations.




