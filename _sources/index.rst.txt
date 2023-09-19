Welcome to Rockable's documentation!
====================================

`Rockable <https://github.com/richefeu/rockable>`_ is a code, written in C++, implementing the classical discrete element method.

The specificities of the code are to handle:

1. Sphero-polyhedral shapes
2. Breakable interfaces

.. note::

   It is designed and developed for academic use (which means that its source code is not publicly downloadable).

.. image:: images/RockableRocks.png
   :width: 192px
   :align: center


Credits
-------

The code was initially developed by *Vincent Richefeu*, at Laboratoire 3SR, to model rockfalls and rock avalanches. This has been done through The PhD work of *Stiven Cuervo* and *Bruna Garcia*, but it actually started before in a code named ``DEMbox`` (no longer maintained).

Then, the breakable interfaces have been implemented during the PhD work of *Marta Stasiak*. A number of improvements have been added at that time thanks to intensive review with *Gael Combe*, Laboratoire 3SR.

New functionalities are being studied thanks to new collaborations of people from CEA, IATE and CNRS. For example,
*Lhassan Amarsid* (CEA) is working on the introduction of periodic boundary conditions, and multi-processor 
computing with domain decomposition. *Farhang Radjai* and students, may introduce new breakable interfaces 
with energy-based criteria. 

Here is the non-exhaustive list of involved people with their main mission: 

- `Vincent Richefeu <vincent.richefeu@univ-grenoble-alpes.fr>`_ (Laboratoire 3SR, UGA): initiated the project
- `Gael Combe <gael.combe@univ-grenoble-alpes.fr>`_ (Laboratoire 3SR, G-INP): mechanical modelling 
- `Lhassan Amarsid <Lhassan.AMARSID@cea.fr>`_ (CEA): mechanical modelling, parallel computing
- `Raphael Prat <Raphael.PRAT@cea.fr>`_ (CEA): parallel computing
- `Jean-Mathieu Vanson <Jean-Mathieu.VANSON@cea.fr>`_ (CEA): mechanical modelling
- `Farhang Radjai <franck.radjai@umontpellier.fr>`_ (LMGC, CNRS): mechanical/physical Mentor
- `Jean-Yves Delenne <jean-yves.delenne@inrae.fr>`_ (IATE, INRAE): mechanical/biological modelling
- `Saeid Nezamabadi <saeid.nezamabadi@umontpellier.fr>`_ (LMGC, UM2): Non Smooth Discrete Element Method (NS-DEM)
- `Patrick Mutabaruka <Patrick.Mutabaruka@ifremer.fr>`_ (Ifremer): coupling with Lattice Boltzmann Method (LBM)


Features
--------

- **Particle Shapes**: the code uses only one 3D shape: sphero-polyhedra or R-shapes. These shapes can be non-convex (with holes if necessary) and have rounded edges and corners (uniform radius per shape).

  .. note:: 
	
	Some other shapes are currently considered for special boundary shapes (sphere, cylinder...) and specifique loadings.
	
- **Boundary Conditions**: any rigid element can be used to apply boundary conditions. It is possible to impose velocity, force, or moment component by component. Some predefined systems with servo-control are also available for complex loading conditions (e.g., loading cycles or controlled pressure).

  .. note:: 
	
	The possibility of applying tri-periodic loading to an assembly is implemented and currently in the testing phase.

- **Parallel Computation**: currently, an OpenMP optimization using compilation flags has been implemented. However, the computational speedup is relatively low. Typically, 8 cores are needed to halve the simulation time (for a dense system with a large number of elements).

- **Documentation**: there is little documentation, although efforts are being made to address this.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickStart
   syntaxConf
   syntaxShapes
   forceLaws
   dissipation
   integrationSchemes
   servoControllers
   preProcessing
