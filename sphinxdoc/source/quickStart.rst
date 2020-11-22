Quick-start guid
================

What is Rockable?
-----------------

Rockable is a DEM code written in C++. The two main specificities of the code are (i) to hold sphero-polyhedral shapes, (ii) to manage breakable interfaces. This technique has been developed in DEMbox which was designed to do more things, with a higher degree of abstraction. Also, it is developed for an academic usage.

Compilation
-----------

The compilation is managed by a makefile that should be similar to that one:

.. code-block:: makefile
   :caption: Makefile
   
   # An option can be disable by adding n before it
   OPTIONS = -DQUAT_ACC -DnFT_CORR -DnROT_MATRIX
   
   UNAME_S := $(shell uname -s)
   ifeq ($(UNAME_S),Darwin)
   # The compiler to be used
   CXX = clang++
   # The list of flags passed to the compiler
   CXXFLAGS = -O3 -Wall -std=c++11 -I ../../common $(OPTIONS)
   # Link flags for OpenGL and glut
   GLLINK = -framework OpenGL -framework GLUT
   else
   # The compiler to be used
   CXX = g++
   # The list of flags passed to the compiler
   CXXFLAGS = -O3 -Wall -std=c++11 -I ../../common $(OPTIONS)
   # Link flags for OpenGL and glut
   GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
   endif
   
   # The list of source files
   SOURCES = Shape.cpp Particle.cpp Interaction.cpp Rockable.cpp DrivingSystem.cpp
   
   # Each cpp file listed below corresponds to an object file
   OBJECTS = $(SOURCES:%.cpp=%.o)
   
   # Each cpp file listed below corresponds to a header file
   HEADERS = $(SOURCES:%.cpp=%.hpp)
   
   # All source files (listed in SOURCES) will be compiled into an object file
   # with the following command
   %.o:%.cpp
     $(CXX) $(CXXFLAGS) -c $<
   
   .PHONY: all clean
   
   all: run see
   
   clean:
     rm -f *.o run see
   
   # The application that runs a simulation
   run: run.cpp $(HEADERS) $(OBJECTS)
     $(CXX) $(CXXFLAGS) -c $<
     $(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)
   
   # The application that visualize the conf files
   see: see.cpp see.hpp $(HEADERS) $(OBJECTS)
     $(CXX) $(CXXFLAGS) -c $< -Wno-deprecated
     $(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLLINK)

As we can see, 3 options can be set at compile time: QUAT_ACC, FT_CORR and ROT_MATRIX. To disable an option it is advised to rename it by adding a letter 'n' at the beginning.

- QUAT_ACC: with this option, the angular position (given by a quaternion) 
  is updated by accounting as prescribed by the chose intergator. 
  Otherwise, the time integration of the quaternions uses a simple explicit Euler scheme

- FT_CORR: with this option, the tangential force is corrected to account for the rotation of the local framework

- ROT_MATRIX: with this option, the computation of the angular accelerations 
  uses rotation matrices instead of quaternion


Depending on the computer you use (Apple or PC with Linux) CXX, CXXFLAGS and GLLINK have to be set in the right section for setting, respectively, the compiler name, the options of compilation, and the options of linkage with openGL.
The rest of the file should normally never be changed.
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

Supposing that the executable named run stands in the same folder as the configuration and shape files, the simulation is launched that way:

.. code-block:: sh
   
   ./run bouncingSphere.txt

If the executable has been compiled with openMP abilities, the number of threads can be set with the option ``-j``, for example:

.. code-block:: sh

   ./run bouncingSphere.txt -j 24

In this particular example, it is clearly not a good idea to use so much threads because the number of particles is to small and the computation duration will be worst.

Visualising the simulations
---------------------------

Normally, the application see has been built as the same time than run. If it is not the case, the compilation can be launched that way:

.. code-block:: sh

   make see

The application ``see`` needs ``freeglut``, the simplest way to use openMP and display 3D things.

