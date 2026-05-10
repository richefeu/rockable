Welcome to Rockable's documentation!
====================================

`Rockable <https://github.com/richefeu/rockable>`_ is a code, written in C++, implementing the classical discrete element method.

The specificities of the code are to handle:

1. Sphero-polyhedral shapes
2. Breakable interfaces

.. note::

   It is designed and developed meanly for academic use.

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
- `Lhassan Amarsid <lhassan.amarsid@cea.fr>`_ (CEA): mechanical modelling, parallel computing
- `Raphael Prat <raphael.prat@cea.fr>`_ (CEA): parallel computing
- `Jean-Mathieu Vanson <jean-mathieu.vanson@cea.fr>`_ (CEA): mechanical modelling
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

- **Parallel Computation**: currently, an OpenMP optimization using compilation flags has been implemented. However, the computational speedup is relatively low. Typically, 8 cores are needed to halve the simulation time (for a dense system with a large number of elements). 

- **HPC Computation**: Rockable can efficiently simulate around 100,000 sphero-polyhedral particles over several weeks. However, simulating millions of particles necessitates high-performance computing (HPC). To maintain the clarity and simplicity of Rockable's source code, HPC capabilities were developed separately. This led to the creation of **ExaDem**, a new HPC code that leverages a subset of Rockable's features for large-scale simulations. The primary objective is to enhance polyhedron-intensive simulations on modern HPC platforms.

  .. note:: 
     The plan is to progressively integrate R-shape functionalities into the **ExaNBody** platform. 
     This will enable the combination of these features with advanced hybrid parallelization techniques, 
     such as MPI + OpenMP and GPU acceleration using CUDA or HIP backends.

- **Periodic boundary conditions**: introducing periodic boundary conditions allows simulations to model systems with repeating patterns or structures, which is often seen in various scientific and engineering contexts. This feature, recently implemented in Rockable enables simulation of systems while reducing edge effects. This is particularly valuable in the context of concurrent double-scale simulations. 

  .. note:: 
     This task is currently handled with Lhassan Amarsid and Duc-Cuong Pham.

- **Non-Smooth Contact Dynamics**: the NSCD approach, long utilized by our research group, is currently being integrated into Rockable. The theoretical modeling is complete, and the next phase involves implementation.

  .. note:: 
     This task is currently handled with Saeid Nezamabadi and Farhang Radjai.

- **Homogeneously compliant particles**: Expanding the framework to handle flexible particles is a significant advancement. This allows simulation of particle materials that can deform under various conditions, such as soft particles, polymers, or biological tissues. The objective is not to determine the kinematic field with a large quantity of degrees of freedom but rather to make certain ones possible.

  .. note:: 
     This task is currently handled in the posdoc of Mukesh Singh Bisht.
   
- **Documentation**: there is little documentation, although efforts are being made to address this.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickStart
   fileOrganisation
   syntaxConf
   syntaxShapes
   preProcessing
   forceLaws
   dissipation
   integrationSchemes
   servoControllers
   visualisation
   tools
   
   
