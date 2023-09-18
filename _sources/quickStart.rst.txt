Quick-start guid
================

What is Rockable?
-----------------

``Rockable`` is a DEM code written in C++11 (tend to become C++17). The two main specificities of the code are (*i*) to hold sphero-polyhedral shapes, (*ii*) to manage breakable interfaces. This technique has been developed in another code named ``DEMbox`` which was designed to do more things, with a higher degree of abstraction. Also, it is designed and developed for an academic usage.

Compilation with cmake
----------------------

The preferred method for building the project is by using CMake, a widely-used build system. Follow these steps to compile the project:

1. Create a new directory called build (you can choose any name) in the project's root directory:

.. code-block:: sh

   mkdir build

2. Navigate into the build directory:

.. code-block:: sh
 
	 cd build
	 
	 
3. Run CMake configuration command from within the build directory. This command sets up the build environment and generates the necessary build files based on the CMake configuration files present in the project:
		
.. code-block:: sh
 
    cmake ..

If you wish to disable the compilation of the ``see`` application, you can pass the ``-DROCKABLE_COMPILE_SEE=OFF`` option to CMake during the configuration step:

.. code-block:: sh

   cmake .. -DROCKABLE_COMPILE_SEE=OFF


To enable the installation with profiling tools, you can use the ``-DROCKABLE_ENABLE_PROFILING=ON`` option during the CMake configuration step:

.. code-block:: sh

   cmake .. -DROCKABLE_ENABLE_PROFILING=ON

These options will set up the CMake build system accordingly to either include or exclude the ``see`` feature and enable or disable the profiling tools during the build process. After configuring CMake with the desired options, you can proceed with the usual build process (*e.g.*, using ``make`` or your preferred build command) to compile the project with the chosen settings.

All the options are listed below:

- ``ROCKABLE_USE_FT_CORR`` (default is ``OFF``): add objectivity correction to tangent forces.
- ``ROCKABLE_ENABLE_PROFILING`` (default is ``ON``): enable time profiling.
- ``ROCKABLE_ENABLE_BOUNDARY`` (default is ``OFF``): enable the special boundaries like Ball or Cylinder.
- ``ROCKABLE_ENABLE_SOFT_PARTICLES``` (default is ``OFF``): enable straining of particles.
- ``ROCKABLE_ENABLE_PERIODIC`` (default is ``OFF``): enable full periodic boundary conditions.
- ``ROCKABLE_COMPILE_SEE`` (default is ``ON``): compile the application to visualize the conf-files.
- ``ROCKABLE_COMPILE_SEE3`` (default is ``OFF``): compile the application to edit graphically the input files.


Compilation with `Makefile`
------------------------------------


The compilation is managed with a classic `Makefile`. Some options can be set at compile time: ``FT_CORR``, and ``COMPONENTWISE_NUM_DAMPING``. To disable an option it is advised to rename it by adding a letter 'n' at the beginning.

- ``FT_CORR``: with this option, the tangential force is corrected to account for the rotation of the local framework
- ``COMPONENTWISE_NUM_DAMPING``: it activates the component-wise numerical damping (so called Cundall damping) 
- ``ROCKABLE_ENABLE_PROFILING``: enable the use of MATools for performance profiling

Depending on the computer you use (Apple or PC with Linux) ``CXX``, ``CXXFLAGS`` and ``GLLINK`` have to be set in the right section for setting, respectively, the compiler name, the options of compilation, and the options of linkage with openGL.
The rest of the file should normally not be changed.
Finally, the applications run and see can be compiled with the following command:

.. code-block:: sh

   make

It is sometimes necessary to remove all object files (.o) together with the compiled applications. this can be made with:

.. code-block:: sh

   make clean


Running a simulation
--------------------


To run a simulation, a configuration file has to be written. The format of such a file is described in the section Syntax for conf-files. We show here a simple example simulating a sphere bouncing on a plan.

.. code-block:: text
   :caption: input.txt
   
   Rockable 20-02-2017
   t 0
   tmax 0.06
   dt 1e-6
   interVerlet 0.01
   interConf 0.01
   
   DVerlet 0.08
   dVerlet 0.02
   density 0 2700
   density 1 2700
   
   forceLaw Avalanches
   knContact 0 1 1e6
   en2Contact 0 1 0.05
   ktContact 0 1 1e7
   muContact 0 1 0.4
   krContact 0 1 1e7
   murContact 0 1 0.0
   
   iconf 0
   nDriven 1
   shapeFile SphereAndPlan.shp
   Particles 2
   Plan 0 0 1 0 -0.05 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
   Sphere 1 0 1 -0.5 0.5 0 3.69 -3.29 0 0 0 0 0.707 0 0.707 0 0 0 -50.52 0 0 0
   
The shape-file as described in the section Syntax for shape-files is a file named SphereAndPlan.shp with the following content:

.. code-block:: text
   :caption: SphereAndPlan.sph
   
   <
   name Plan
   radius 0.05
   preCompDone y
   nv 4
   2 0 0.5
   2 0 -0.5
   -2 0 -0.5
   -2 0 0.5
   ne 4
   0 1
   1 2
   2 3
   3 0
   nf 1
   4 0 1 2 3
   obb.extent 2.0 0.05 0.5
   obb.e1 1 0 0
   obb.e2 0 1 0
   obb.e3 0 0 1
   obb.center 0 0 0
   volume 1
   I/m 1 1 1
   >
   
   <
   name Sphere
   radius 0.08
   preCompDone y
   nv 1
   0 0 0
   ne 0
   nf 0
   obb.extent 1 1 1
   obb.e1 1 0 0
   obb.e2 0 1 0
   obb.e3 0 0 1
   obb.center 0 0 0
   volume 0.004021
   I/m 0.00493333 0.00493333 0.0032
   >

Supposing that the executable named ``rockable`` stands in the same folder as the configuration and shape files, the simulation is launched that way:

.. code-block:: sh
   
   ./rockable bouncingSphere.txt

If the executable has been compiled with openMP abilities, the number of threads can be set with the option ``-j``, for example:

.. code-block:: sh

   ./rockable bouncingSphere.txt -j 24

In this particular example, it is clearly not a good idea to use so much threads because the number of particles is to small and the computation duration will be worst.

The verbosity of logs is set with a number that way:

.. code-block:: sh

  ./rockable bouncingSphere.txt -v 6

Highest number corresponds highest verbosity: ``trace`` = 6, ``debug`` = 5, ``info`` = 4, ``warn`` = 3, ``err`` = 2, ``critical`` = 1, ``off`` = 0

If the files produced by a computation (``conf*``, ``kineticEnergy.txt``, ``perf.txt``, and ``staticBalance.txt``) have to be deleted, ``rockable`` can do the job.

.. code-block:: sh

  ./rockable -c


Visualising the simulations
---------------------------

Normally, the application ``see`` has been built as the same time than ``rockable``. 
The application ``see`` needs ``freeglut``, the simplest way to use openMP and display 3D things.

