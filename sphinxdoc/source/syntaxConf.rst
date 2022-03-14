
Format of configuration files (``conf-file``)
=============================================

The conf-files are the files that hold the whole configuration at a given time. They are used:

 1. for defining the initial configuration and parameters of a simulation, 
 2. for running some preprocessing commands,
 3. for saving the history of a simulation. The keywords are defined in the following.

Header
------

A conf-file always starts with the header: ``Rockable dd-mm-yyyy`` (*e.g.*, ``Rockable 20-02-2017``). 
Each time a noticeable change is made in the format, the date of this change is also changed in the header of the file. It is used as the version of the format.

Timing
------

- ``t`` (*double*) **value**  
  Current time

- ``tmax`` (*double*) **value**  
  Maximum time (the simulation will end when time is ``tmax``)

- ``dt`` (*double*) **value**  
  Time step increments


Neighbor list (``NL``)
----------------------

- ``interVerlet`` (*double*) **value**  
  Elapsed time between each rebuilding of the neighbor list

- ``DVerlet`` (*double*) **value**  
  Distance to define if two sphero-polyhedra are neighbors

- ``dVerlet`` (*double*) **value**  
   Distance to define if two sub-elements (sphere for vertices, tubes for edges, 
   or thick 3D polygons for faces), between two sphero-polyhedra, are neighbors

- ``dynamicUpdateNL`` 0/1  
  If the dynamic update of the neighbor list is activated (set to 1), 
  then an update will be made if the maximum distance of a body since the last update becomes 
  larger than ``dispUpdateNL``. 
  An update will also be made when the maximum rotation becomes larger than ``angleUpdateNL``.

   - ``dispUpdateNL`` (*double*) **distance**

   - ``angleUpdateNL`` (*double*) **angleDegree**

Configuration backups
---------------------

- ``interConf`` (*double*) **value**  
  Elapsed time between each backup of the configuration

- ``iconf`` (*double*) **value**  
  Number of the current configuration

  
Options
-------

- ``AddOrRemoveInteractions`` (*string*) **Option** where **Option** can be ``bruteForce`` or ``OBBtree``. Default is ``bruteForce``

- ``UpdateNLStrategy`` (*string*) **Option** where **Option** can be ``bruteForce`` or ``linkCells``. Default is ``bruteForce``

  In case of ``linkCells`` strategy, the vector ``cellMinSizes`` needs to be set. Each value is the minimum size of a cell in the corresponding direction. In addition, ``boxForLinkCellsOpt`` is a flag (0 or 1) to state if the first driven bodies are part of the overall bounding box (that will be split in to cells)
  
  -  ``cellMinSizes`` (*double*) **xmin** (*double*) **ymin** (*double*) **zmin**

  -  ``boxForLinkCellsOpt`` (0 | 1)

- ``Integrator`` (string) **Option** where **Option** can be ``Euler``, ``velocityVerlet``, ``Beeman`` or ``RungeKutta4``. Default is ``velocityVerlet`` 

Library of particle shapes
--------------------------

- ``shapeFile`` (*string*) **path**  
  Path of the file that define the shapes used. 
  The format to define a shape is explained in a separated document (``SyntaxShape.pdf``)

Particles
---------

- ``density`` (*int*) **groupNumber** (*double*) **density**

- ``Particles`` (*int*) **numberOfParticles**  

  Then these entries are repeated: (*string*)shapeName (*int*) **group** (*int*) **cluster** (*double*) **homothety** (*vec3r*) **position** (*vec3r*) **velocity** (*vec3r*) **acceleration** (*quat*) **angularPosition** (*vec3r*) **angularVelocity** (*vec3r*) **angularAcceleration**


Interactions
------------

- ``Interactions`` (*int*) **numberOfInteractions**  

  Then these **numberOfInteractions** entries are repeated: 
  (*int*) **i** (*int*) **j** (*int*) **type** (*int*) **isub** (*int*) **jsub** (*vec3r*) **n**  (*double*) **dn**
  (*vec3r*) **position** (*vec3r*)* *relativeVelocity** (*double*) **fn** (*vec3r*) **ft** (*vec3r*) **mom**
  (*double*) **viscousDampingValue**  


    **type** is either 0 for vertex-vertex, 1 for vertex-edge, 2 for vertex-face, or 3 for edge-edge
    **relativeVelocity** is the velocity of body **j** relative to body **i** at the contact point
    **n** is oriented from **j** to **i**

Force laws
----------

- ``forceLaw`` (*string*) **Name**  
  Select a model for the computation of forces. For possible **Name** see :ref:`Force-laws`

Time-integration scheme
-----------------------

- ``xxxxx`` (*string*) **Name**  
  Select a scheme for the time integration. For possible **Name** see :ref:`IntegrationSchemes`

Dissipation
-----------

There are several dissipation strategies.

Loading
-------

- ``nDriven`` (*int*) **value**
  **value** is the number of bodies, at the beginning of the list, that are not free to move. 
  By default, the **nDriven** first bodies are fixed (all velocities imposed to zero), 
  but if we want to set a velocity or a force/moment, some commands have to be added 
  in a file named ``drivingSystem.txt``.

File drivingSystem.txt
----------------------

- ``Control`` (*string*) **mode** (*int*) **bodyNumber** (*double*) **value**
  where **mode** is either ``_x_Vel_``, ``_y_Vel_``, ``_z_Vel_``, ``_xrot_Vel_``, ``_yrot_Vel_``, 
  ``_zrot_Vel_``, ``_x_For_``, ``_y_For_``, ``_z_For_``, ``_xrot_Mom_``, ``_yrot_Mom_``, or ``_zrot_Mom_``  
  
  .. note:: There are 2 more mode keywords for which the single **value** has to be replaced by 3 values
            (a vector of 3 components):  ``_xyzrot_Vel_`` and ``_xyzrot_Mom_``
  
- ``Servo`` (*string*) **servoName** <*PARAMETERS*>

Processing commands
-------------------

* ``stickVerticesInClusters`` (*double*) **Epsilon** 
  This command will add glued interfaces between bodies having the same cluster identifier. 
  Only bonds between vertices (spheres) are added when the distance is less than **Epsilon**.

* ``stickClusters`` (*double*) **Epsilon**  
  This command will add glued interfaces between bodies having different cluster identifier. 
  Bonds are added when the distance is less than **Epsilon**.
  
* ``setAllVelocities`` (*vec3r*) **velocity**
  Set the velocity vector of all particles (that are not driven) to the prescribed value **velocity**.

Data Extractors
---------------

- ``DataExtractor`` (*string*) **ExtractorName** <*PARAMETERS*>  
  The list of <*PARAMETERS*> depends on the ExtractorName. They are listed below.



