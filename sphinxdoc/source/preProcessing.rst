.. _prePro:

Pre-processing commands
=======================

Pre-processing commands are keywords of the conf-file that **act on the system
once it has been read**, rather than describing it. They are the way to glue
bodies together, to randomise velocities, or to duplicate a group of particles,
without having to generate all of it explicitly.

.. important::

   They act on a system that is already set, so they belong at the **end** of
   the input file, after ``Particles`` and, when relevant, after
   ``Interactions``. A command placed before the particles are read would find
   nothing to work on.

.. note::

   A pre-processing command is executed **once**, when the file is read, and it
   is never written back into the conf-files that the computation saves. The
   dumps hold the *result* of the command: the interfaces it created, or the
   velocities it set. Restarting from a dump therefore does not re-run it,
   which is what you want.

.. contents::
   :local:
   :depth: 2


Creating glued interfaces
-------------------------

These commands build the ``Interfaces`` of the sample. Two bodies are glued when
they are closer than a distance **Epsilon**, and the resulting bonds are what a
breakable force law (``StickedLinks``, ``BCM``) acts on.

A bond is said to be **inner** when the two bodies share the same ``cluster``
number, and **outer** otherwise. This is what selects the ``*InnerBond`` or
``*OuterBond`` parameters.

``stickVerticesInClusters`` (*double*) **Epsilon**
    Glue the bodies that have the **same** cluster identifier. Only
    vertex-vertex bonds (sphere to sphere) are created, when the distance is
    less than **Epsilon**.

``stickVerticesInClustersMoments`` (*double*) **Epsilon**
    Same as above, with the transmission of moments enabled at the bonds. Use
    this one when the assembly has to resist bending rather than only traction.

``stickClusters`` (*double*) **Epsilon**
    Glue bodies belonging to **different** clusters. All the bond types are
    created, not only vertex-vertex, when the distance is less than
    **Epsilon**.

``stickBCM`` (*double*) **Epsilon**
    Build the interfaces required by the ``BCM`` force law. Beyond creating the
    bonds, it computes the **area** of each interface, which is what the
    energy-based rupture criterion needs, and copies ``knInnerBond``,
    ``ktInnerBond`` and ``gcInnerBond`` into each interface. All the interfaces
    it creates are flagged as inner.

.. tip::

   ``glue_with_walls yes``, placed before these commands, extends them to the
   driven bodies, so that a sample can be glued to its container.

.. note::

   Choosing **Epsilon** is a geometric matter: it must be large enough to catch
   the neighbours that should be bonded, and small enough not to bond bodies
   that merely pass close to each other. A value of the order of the Minkowski
   radius is a reasonable starting point.


Tuning the interfaces
---------------------

These commands modify interfaces that already exist, so they must come **after**
one of the sticking commands.

``copyParamsToInterfaces`` (*string*) **inner|outer**
    Copy the parameters of the group-pair table into each interface, so that
    every interface then carries its own copy. Only the interfaces of the
    selected kind are touched: ``inner`` visits the interfaces whose two bodies
    share a cluster, ``outer`` the others.

    This command sets ``ParamsInInterfaces`` to 1 by itself, so the copied
    values are saved in the conf-files and read back on a restart.

    .. warning::

       The two branches do not copy the same set. ``outer`` copies
       :math:`k_n, k_t, k_r, f_n^0, f_t^0, M_0` and the exponent, whereas
       ``inner`` copies only :math:`k_n, k_t, f_n^0, f_t^0` and the exponent:
       the rolling stiffness and the moment threshold are **not** copied for
       inner bonds, and keep whatever the interface already held. Set them
       explicitly if your model relies on them.

``setStiffnessRatioInterfaces`` (*double*) **ratio**
    Set :math:`k_t = \text{ratio} \times k_n` in every interface.

``setVariableStickParams`` (*string*) **paramName** (*string*) **inner|outer** (*double*) **lambda** (*int*) **m** (*int*) **timeSeeded**
    Draw one interface parameter from a **Weibull** distribution of scale
    ``lambda`` and modulus ``m``, independently for each interface. This is the
    usual way of introducing a controlled scatter of strength, and hence a
    progressive rather than simultaneous failure.

    ``timeSeeded`` set to 1 seeds the random generator on the clock, so that two
    runs differ; set to 0, the same sample is reproduced exactly.

.. important::

   These commands write into the interfaces themselves, and the force laws only
   read what an interface carries when ``ParamsInInterfaces`` is 1.
   ``copyParamsToInterfaces`` sets that flag on its own; the two others do not,
   so set ``ParamsInInterfaces 1`` in the conf-file if you use them alone.


Setting the kinematics
----------------------

The velocities of the **driven** bodies are never modified by these commands.

``setAllVelocities`` (*vec3r*) **velocity**
    Set the velocity vector of all the free particles to the prescribed vector.

``randomlyOrientedVelocities`` (*double*) **velocityMagnitude**
    Give every free particle a velocity of the prescribed magnitude, in a
    direction drawn uniformly at random. Each particle gets its own direction.

``randomlyOrientedVelocitiesClusters`` (*double*) **velocityMagnitude** (*int*) **opt**
    Same, but a single direction is drawn per **cluster**, so that the bodies of
    a cluster keep moving together instead of being torn apart at the first
    step.

    ``opt`` set to 1 forces the vertical component downwards
    (:math:`v_y \leftarrow -|v_y|`), which is what a release of blocks above a
    slope needs; set to 0, the direction stays uniform over the sphere.


Populating the sample
---------------------

``homothetyRange`` (*int*) **ifirst** (*int*) **ilast** (*double*) **hmin** (*double*) **hmax** (*int*) **timeSeeded**
    Give the particles from ``ifirst`` to ``ilast`` a homothety drawn uniformly
    in :math:`[h_\text{min}, h_\text{max}]`. This is how a size dispersion is
    introduced without defining one shape per size.

    ``timeSeeded`` behaves as in ``setVariableStickParams``.

    .. note::

       The command recomputes the mass and the inertia of each particle it
       touches, from the new homothety, the volume of its shape and the density
       of its group. The ``density`` of every group concerned must therefore
       already be defined when the command is executed, which is the case as
       long as it sits at the end of the file.

``particlesClonage`` (*int*) **ifirst** (*int*) **ilast** (*vec3r*) **translation**
    Duplicate the particles from ``ifirst`` to ``ilast``, translated by the
    given vector. Useful to build a periodic-looking pattern, or to stack
    several copies of a prepared block.

    The clones are appended at the end of the particle list, and they are given
    **new cluster numbers**, continuing after the highest one in use. A cloned
    block is therefore a separate cluster: bonds inside a clone stay inner,
    while the clone is not glued to its original.

    .. note::

       The clones are appended after the file has been read, so the particle
       count of the input file no longer matches the number of bodies in
       memory. The next conf-file written by the computation records the new
       count, and is a consistent input on its own.


A typical bonded sample
-----------------------

.. code-block:: text
   :caption: input.txt (end of the file)

   ...
   ParamsInInterfaces 1
   Particles 512
   ...

   stickVerticesInClusters 0.05
   copyParamsToInterfaces inner
   setVariableStickParams fn0InnerBond inner 1.0e4 5 0

The sample is glued cluster by cluster, each interface then receives its own
copy of the parameters, and the tensile strength of the bonds is finally
scattered with a Weibull law of modulus 5.
