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

- ``UpdateNLStrategy`` (*string*) **Option**

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

- ``Control`` (*string*) **mode** (*int*) **bodyNumber** (*double*) **value**

  Use **mode** to set a velocity or force/moment for a specific body.
  
  .. warning:: 
	
	There are two additional mode keywords for which the single **value** has to be replaced by three values 
	(a vector of three components): ``_xyzrot_Vel_`` and ``_xyzrot_Mom_``.

- ``Servo`` (*string*) **servoName** <*PARAMETERS*>

  Set the parameter list depending on the selected servo (see :ref:`Servo-controllers`).


Pre-processing Commands
-----------------------

``Rockable`` provides several commands for performing preprocessing tasks.
These commands are typically entered at the end of an input ``conf-file``, after the definition of particles
and interactions. For details on these commands, refer to: :ref:`prePro`


Data Extractors
---------------

- ``DataExtractor`` (*string*) **ExtractorName** <*PARAMETERS*>

  The parameter list depends on the **ExtractorName**. Details are listed below.


