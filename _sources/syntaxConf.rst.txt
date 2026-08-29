.. _syntaxConf:

Format of configuration files (conf-file)
==========================================

The conf-files are the files that hold the whole configuration at a given time. They are used for the following purposes:

1. Defining the initial configuration and parameters of a simulation.
2. Running some preprocessing commands.
3. Saving periodically the history of a simulation.


Header
------

A conf-file always starts with the following header format: ``Rockable dd-mm-yyyy`` (e.g., ``Rockable 29-11-2018``). 
The header includes the date (``dd-mm-yyyy``) corresponding to the version of the format. 
Whenever a noticeable change is made in the format, the date in the header is updated to reflect the new version. 
The version date is defined in the preprocessor define of ``CONF_VERSION_DATE`` in the source code.


Keywords
--------

The conf-files contain specific keywords that define various aspects of the configuration. 
These keywords are used throughout the file to specify different settings and parameters.


Usage
-----

These conf-files are essential for setting up and managing simulations. They are used for:

1. Initializing the simulation with specific configurations and parameters.
2. Executing preprocessing commands to prepare the simulation environment.
3. Periodically saving simulation history to track its progress over time.

The format of the conf-files allows users to define and control various aspects of the simulation accurately.

.. tip:: 

   It is crucial to maintain consistency in the format of conf-files and keep track of changes made 
   to the header date whenever modifications are introduced. This practice ensures clarity and compatibility 
   between different versions of the configuration files and the codebase.


Timing
------

- ``t`` (*double*) **value**

  Current time.

- ``tmax`` (*double*) **value**

  Maximum time. The simulation will end when time reaches ``tmax``.

- ``dt`` (*double*) **value**

  Time step increments.

Neighbor List (NL)
------------------

- ``interVerlet`` (*double*) **value**

  Elapsed time between each rebuild of the neighbor list.

- ``DVerlet`` (*double*) **value**

  Distance used to define if two sphero-polyhedra are neighbors. 
  This length is added to the Oriented Bounding Boxes (OBBs) before testing for overlap. 
  In other words, half of this length is added to each side of the OBBs.

- ``dVerlet`` (*double*) **value**

  Distance used to define if two sub-elements (sphere for vertices, tubes for edges, 
  or thick 3D polygons for faces) between two sphero-polyhedra are neighbors.

- ``dynamicUpdateNL`` (0 | 1)

  If dynamic update of the neighbor list is activated (set to 1), 
  the list will be updated if the maximum distance a body has moved since the last update 
  becomes larger than ``dispUpdateNL``, or if the maximum rotation becomes larger than ``angleUpdateNL``. 
  These dynamic updates will not affect the regular updates (every ``interVerlet``).

  Additional parameters to be set:

  - ``dispUpdateNL`` (*double*) **distance**
  - ``angleUpdateNL`` (*double*) **angleDegree**


Parallel mode
-------------

The ``parallel_mode`` parameter controls the parallel execution strategy used by Rockable for force computation and interaction handling.
It defines how work is distributed across threads and how memory conflicts are handled during updates.

Available options
~~~~~~~~~~~~~~~~~

``DefaultParallelMode``
    Uses the default execution strategy defined by Rockable.  
    Typically mapped to the most stable or recommended mode depending on the build configuration.

``InteractionBuffer``
    Uses an intermediate buffer to store interactions before applying updates.  

``CellMutex``
    Assigns one thread per cell and protects shared data using mutexes.  
    Simple and safe, but may introduce synchronization overhead.

``WaveMethod``
    Uses a wave-based decomposition with 27 waves.  
    One thread processes one cell without mutexes, relying on a structured traversal order to avoid conflicts.

``WaveMethodBlock``
    Similar to ``WaveMethod`` but operates on blocks of at least 8 cells.  
    Uses an 8-wave decomposition and improves cache efficiency while maintaining lock-free execution.

Input keyword
~~~~~~~~~~~~~

The input parameter ``parallel_mode`` is parsed from string values.

Accepted keywords are **exactly identical** to the enum names:

- ``DefaultParallelMode``
- ``InteractionBuffer``
- ``CellMutex``
- ``WaveMethod``
- ``WaveMethodBlock``

Example
~~~~~~~

.. code-block:: yaml

    parallel_mode WaveMethodBlock


Configuration Backups
---------------------

As mentioned earlier, the "conf-files" store a configuration (i.e., the entire dataset of a simulation state). 
The format described here can be used as both input and output of a simulation.

- ``interConf`` (*double*) **value**

  Elapsed time between each backup of the configuration.

- ``iconf`` (*double*) **value**

  Number of the current configuration. This number is used to name the conf-file.

	
Computation options
-------------------


- ``AddOrRemoveInteractions`` (*string*) **Option**

  Choose the method for handling interactions between particles (sphero-polyhedra).

  - **Options**: ``bruteForce`` (default) or ``OBBtree``
  - **Description**: The interactions between particles involve different types of interaction, such as sphere-sphere, sphere-tube, sphere-polygon, and tube-tube. The best strategy to be used depends on the complexity of the involved shapes.

- ``UpdateNL`` (*string*) **Option**

  Choose the strategy for updating neighbor lists.

  - **Options**: ``bruteForce`` (default) or ``linkCells``
  - **Description**: In the case of the ``linkCells`` strategy, the following additional settings need to be configured:

    - ``cellMinSizes`` (*double*) **xmin** (*double*) **ymin** (*double*) **zmin**: Set the minimum size of a cell in each direction.
    - ``boxForLinkCellsOpt`` (0 | 1): A flag to determine if the first driven bodies are part of the overall bounding box, which will be split into cells.

- ``Integrator`` (*string*) **Option**

  Choose the time-integration scheme to be used.

  - **Options**: ``Euler``, ``velocityVerlet`` (default), ``Beeman``, or ``RungeKutta4``
  - **Description**: This option determines the time-integration method for the simulation. For more details, 
	see :ref:`IntegrationSchemes`.


Library of Particle Shapes
--------------------------

- ``shapeFile`` (*string*) **path**

  Path of the file that defines the shapes used.
  The format to define a shape is explained here: :ref:`syntaxShape`


Particles
---------

- ``density`` (*int*) **groupNumber** (*double*) **density**

  Set the density (in kilograms per cubic meter) for particles belonging to a given group number.

- ``Particles`` (*int*) **numberOfParticles**

  The following entries are repeated for each particle:
  (*string*)shapeName (*int*) **group** (*int*) **cluster** (*double*) **homothety** (*vec3r*) **position** 
  (*vec3r*) **velocity** (*vec3r*) **acceleration** (*quat*) **angularPosition** (*vec3r*) **angularVelocity** 
  (*vec3r*) **angularAcceleration**


Interactions
------------

- ``Interactions`` (*int*) **numberOfInteractions**

  The following **numberOfInteractions** entries are repeated for each interaction:
  (*int*) **i** (*int*) **j** (*int*) **type** (*int*) **isub** (*int*) **jsub** (*vec3r*) **n** (*double*) **dn**
  (*vec3r*) **position** (*vec3r*) **relativeVelocity** (*double*) **fn** (*vec3r*) **ft** (*vec3r*) **mom**
  (*double*) **viscousDampingValue**

  - **type**: 0 for vertex-vertex, 1 for vertex-edge, 2 for vertex-face, or 3 for edge-edge.
  - **relativeVelocity**: The velocity of body **j** relative to body **i** at the contact point.
  - **n**: Vector oriented from **j** to **i**.


Force Laws
----------

- ``forceLaw`` (*string*) **Name**

  Select a model for the computation of forces. For possible **Name**, see :ref:`Force-laws`.


Time-Integration Scheme
-----------------------

- ``Integrator`` (*string*) **Name**

  Select a scheme for time integration. For possible **Name**, see :ref:`IntegrationSchemes`.


Dissipation
-----------

There are several dissipation strategies that can be used (see :ref:`Dissipation`).


Loading
-------

- ``nDriven`` (*int*) **Value**

  Set the number of bodies, at the beginning of the list, that are not free to move. 
  By default, the **nDriven** first bodies are fixed (all velocities imposed to zero), 
  but if we want to set a velocity or a force/moment, some commands have to be added 
  in a file named ``drivingSystem.txt``.


File drivingSystem.txt
----------------------

The bodies driven with ``nDriven`` are piloted from a separate file named
``drivingSystem.txt``, which holds the ``Control`` and ``Servo`` entries.
It has a page of its own: :ref:`drivingSystem`, and the servos are listed in
:ref:`Servo-controllers`.


Pre-processing Commands
-----------------------

``Rockable`` provides several commands for performing preprocessing tasks.
These commands are typically entered at the end of an input ``conf-file``, after the definition of particles
and interactions. For details on these commands, refer to: :ref:`prePro`


Data Extractors
---------------

Measurements evaluated during the computation are declared in a file named
``dataExtractors.txt``; see :ref:`dataExtractors`. The keyword
``DataExtractor`` is also accepted inside a conf-file, but it is deprecated
because it is not written back into the saved conf-files.


Body Forces
-----------

- ``BodyForce`` (*string*) **Name** <*PARAMETERS*>

  Apply a force computed from the state of each body alone. Only one body force
  can be active at a time. For possible **Name** and their parameters, see
  :ref:`bodyForces`.


Body properties and environment
-------------------------------

- ``density`` (*int*) **groupNumber** (*double*) **value**

  Density, in kg/m³, of the bodies of a given group. The mass and the inertia of
  a particle are computed from this value, the volume of its shape and its
  homothety, at the moment the particle is read. ``density`` must therefore come
  **before** ``Particles``.

- ``gravity`` (*vec3r*) **vector**

  Gravity acceleration vector. It is a vector, not a magnitude, so an inclined
  gravity is the simplest way to tilt a whole sample without moving anything.

- ``precision`` (*int*) **nbDigits**

  Number of significant digits used when writing the conf-files. Increasing it
  makes a restart closer to a continuous run, at the cost of larger files.


Interfaces and bonds
--------------------

- ``ParamsInInterfaces`` (*int*) **0|1**

  When set to 1, each interface carries its own copy of the bond parameters,
  written on its line in the conf-file, instead of taking them from the
  group-pair table. This is required by the pre-processing commands that
  scatter the bond properties (see :ref:`prePro`).

- ``Interfaces`` (*int*) **numberOfInterfaces**

  Opens the list of glued interfaces. Each entry holds, on one line:

  ``<i> <j> <nbBonds> <dn0>``

  then, only when ``ParamsInInterfaces`` is 1,

  ``<kn> <kt> <kr> <fn0> <ft0> <mom0> <power> <Gc>``

  and finally ``nbBonds`` triplets ``<type> <isub> <jsub>`` identifying the
  bonded sub-interactions among those already read in ``Interactions``.

  .. warning::

     A bond that cannot be matched with an existing interaction produces a
     warning and the **whole interface** is dropped. This is why ``Interactions``
     must be read before ``Interfaces``.

- ``glue_with_walls`` (*string*) **yes|no**

  When set to yes, the sticking pre-processing commands also glue the free
  bodies to the driven ones. Accepted affirmative forms are ``yes``, ``YES``,
  ``y``, ``Y`` and ``1``.

- ``initSpringJoint`` (*int*) **ibody** (*vec3r*) **ipos0** (*int*) **jbody** (*vec3r*) **jpos0** (*double*) **stiffness**

  Add a linear spring joining two bodies. The anchor points ``ipos0`` and
  ``jpos0`` are given in the **body frames**, so they follow the bodies as they
  rotate. Unlike an interface, a spring joint never breaks.


Interaction parameters
----------------------

All of them follow the same pattern, described in :ref:`Force-laws`:

``<parameter> <group1> <group2> <value>``

The table is symmetric: setting ``0 2`` also sets ``2 0``. A pair that is never
defined keeps a value of zero, so **every pair of groups that can meet must be
given the parameters of the selected force law**.

.. list-table:: Contact
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Meaning
   * - ``knContact``
     - Normal stiffness :math:`k_n`.
   * - ``ktContact``
     - Tangential stiffness :math:`k_t`.
   * - ``muContact``
     - Coulomb friction coefficient :math:`\mu`.
   * - ``krContact``
     - Rolling stiffness :math:`k_r`. The resistant moment is computed only when
       it is strictly positive.
   * - ``murContact``
     - Rolling resistance coefficient :math:`\mu_r`, used as a **length**.
   * - ``en2Contact``
     - Squared normal restitution coefficient :math:`e_n^2`.
   * - ``en2ContactFromViscRate``
     - The same quantity, given the other way round: the value is the viscous
       damping rate :math:`\alpha_n`, and :math:`e_n^2` is deduced from it.
   * - ``viscTContact``
     - Tangential viscosity, read only by the ``GeoVisc`` law.

.. list-table:: Bonds
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Meaning
   * - ``knInnerBond``, ``knOuterBond``
     - Normal stiffness of the bond.
   * - ``ktInnerBond``, ``ktOuterBond``
     - Tangential stiffness of the bond.
   * - ``krInnerBond``, ``krOuterBond``
     - Rolling stiffness of the bond.
   * - ``en2InnerBond``, ``en2OuterBond``
     - Squared normal restitution of the bond.
   * - ``fn0InnerBond``, ``fn0OuterBond``
     - Normal (tensile) strength :math:`f_n^0`.
   * - ``ft0InnerBond``, ``ft0OuterBond``
     - Tangential strength :math:`f_t^0`.
   * - ``mom0InnerBond``, ``mom0OuterBond``
     - Moment strength :math:`M_0`.
   * - ``powInnerBond``, ``powOuterBond``
     - Exponent :math:`p` of the rupture criterion.
   * - ``gcInnerBond``, ``gcOuterBond``
     - Fracture energy :math:`G_c`, read by the ``BCM`` law.

A bond is **inner** when the two bodies share the same ``cluster`` number, and
**outer** otherwise.

- ``ContactPartnership`` (*string*) **model**

  Distributes the interaction parameters over the sub-contacts of a same pair of
  bodies, so that two bodies touching through many sub-contacts are not stiffer
  than two bodies touching through one. **model** is ``None`` (default),
  ``NumberWeight``, ``OverlapWeight`` or ``SurfaceWeight``.

- ``preventCrossingLength`` (*double*) **value**

  When strictly positive, stiffens the normal repulsion at large overlap, to
  prevent bodies from crossing each other. It is written in the conf-files only
  when it is non-zero.


Time-dependent parameters
-------------------------

The ``Tempo`` keyword plugs a value onto a function of time, re-evaluated at
every step. Two profiles exist:

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - Profile
     - Behaviour
   * - ``Range <t1> <t2> <v1> <v2>``
     - ``v1`` while :math:`t \in [t_1, t_2]`, ``v2`` outside.
   * - ``Ramp <t1> <t2> <v1> <v2>``
     - ``v1`` before :math:`t_1`, linear interpolation up to ``v2`` at
       :math:`t_2`, then ``v2``.

Three targets can be driven:

- ``Tempo NDCoeff`` <*PROFILE*>

  Drives the numerical damping coefficient.

- ``Tempo Inter`` (*string*) **paramName** (*int*) **g1** (*int*) **g2** <*PROFILE*>

  Drives one interaction parameter of a group pair. Both ``g1 g2`` and
  ``g2 g1`` are updated.

- ``Tempo Body`` (*string*) **propName** (*int*) **group** <*PROFILE*>

  Drives one body property, such as ``density``, of a group.

.. code-block:: text
   :caption: input.txt

   # relax the friction progressively over the first half second
   Tempo Inter muContact 0 1 Ramp 0.0 0.5 0.9 0.2

   # damp strongly, but only at the beginning
   Tempo NDCoeff Range 0.0 0.1 0.05 0.0


Velocity barriers
-----------------

An alternative to numerical damping, described in :ref:`Dissipation`.

- ``VelocityBarrier`` (*double*) **value** and ``VelocityBarrierExponent`` (*double*) **value**

  Cap on the translational velocities, and the exponent of the barrier function.

- ``AngularVelocityBarrier`` (*double*) **value** and ``AngularVelocityBarrierExponent`` (*double*) **value**

  The same, for the rotational velocities.


Optional features
-----------------

These keywords are only understood when Rockable has been compiled with the
corresponding CMake option. Otherwise they are accepted but print a message
saying the feature was not enabled; see :ref:`commandLine`.

- Periodic cell: ``usePeriodicCell``, ``h``, ``vh``, ``ah``, ``mh``, ``dh``,
  ``cellVelocityCorrection``, ``cellMomentumCorrection``, ``useKineticStress``.
  See :ref:`periodicBoundaryConditions`.

- Deformable particles: ``useSoftParticles``. See :ref:`softParticles`.


Deprecated keywords
-------------------

- ``separator`` (*string*) **keyword**

  Sets the separator used when writing the conf-files. It is marked for removal
  in the source and should not be used.

- ``OBBtreeLevel``, in the shape files, is likewise deprecated.
