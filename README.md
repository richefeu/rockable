<p align="center">
<img src="./sphinxdoc/source/images/RockableRocks.png" width="50%"/>
</p>

## What is `Rockable`?  

Rockable is a DEM code written in C++ by <vincent.richefeu@3sr-grenoble.fr>. The two main specificities of the code are (_i_) to hold sphero-polyhedral shapes, (_ii_) to manage breakable interfaces. It is developed for an **academic usage**. This means that the code is not intended to be a tool for all purposes. It can easily be used to do what it is designed for, but to extend it, it is necessary to master both the model (DEM, complex shapes and interaction laws) and its implementation (data structure). The benefit of a good understanding is to avoid a "hacking" that would eventually limit the developed possibilities. In other words, the design of the code (neither too specific nor too general) is intended to avoid any tendency towards a single thought.

The use of the code is not interfaced by any tool (like lua, python or any graphical interface) to facilitate its use, except the input format as described in the documentation. This makes it particularly streamlined and greatly facilitates its integration with other calculation codes. It is in this sense that Rockable is qualified of "academic code".

### Source tree 

* `doc`: place to generate the doxygen html documentation of the source files of Rockable
* `sphinxdoc`: user documentation (sphinx with ReStructuredText)
* `examples`: examples for usage tutorials or for testing features
* `prepro`: some pre-processing tools
* `src`: C++ source files
* `test`: regression test files

### Folders created at building stage

* `deps`: source files for Rockable dependencies fetched by `cmake`
* `BUILD`: compilation files of the code (created by the script `install_rockable.sh`)
* `INSTALL:` binaries of Rockable routines (created by the script `install_rockable.sh`)

## Credits

The code was initially developed by *Vincent Richefeu*, at Laboratoire 3SR, to model rockfalls and rock avalanches. This has been done through The PhD work of *Stiven Cuervo* and *Bruna Garcia*, but it actually started before in a code named ``DEMbox`` (no longer maintained).

Then, the breakable interfaces have been implemented during the PhD work of *Marta Stasiak*. A number of improvements have been added at that time thanks to intensive review with *Gael Combe*, Laboratoire 3SR.

New functionalities are being studied thanks to new collaborations of people from CEA, IATE and CNRS. For example,
*Lhassan Amarsid* (CEA) is working on the introduction of periodic boundary conditions, and multi-processor 
computing with domain decomposition. *Farhang Radjai* and students, may introduce new breakable interfaces 
with energy-based criteria. 

Here is the non-exhaustive list of involved people with their main mission: 

* **Vincent Richefeu** <vincent.richefeu@univ-grenoble-alpes.fr> (Laboratoire 3SR, UGA): initiated the project
* **Gael Combe** <gael.combe@univ-grenoble-alpes.fr> (Laboratoire 3SR, G-INP): mechanical modelling 
* **Lhassan Amarsid** <Lhassan.AMARSID@cea.fr> (CEA): mechanical modelling, parallel computing
* **Raphael Prat** <Raphael.PRAT@cea.fr> (CEA): parallel computing
* **Jean-Mathieu Vanson** <Jean-Mathieu.VANSON@cea.fr> (CEA): mechanical modelling
* **Farhang Radjai** <franck.radjai@umontpellier.fr> (LMGC, CNRS): mechanical/physical Mentor
* **Jean-Yves Delenne** <jean-yves.delenne@inrae.fr> (IATE, INRAE): mechanical/biological modelling
* **Saeid Nezamabadi** <saeid.nezamabadi@umontpellier.fr> (LMGC, UM2): Non Smooth Discrete Element Method (NS-DEM)
* **Patrick Mutabaruka** <Patrick.Mutabaruka@ifremer.fr> (Ifremer): coupling with Lattice Boltzmann Method (LBM)


## Features

* Particle Shapes: the code uses only one 3D shape: sphero-polyhedra or R-shapes. These shapes can be non-convex (with holes if necessary) and have rounded edges and corners (uniform radius per shape).

  > [!NOTE]
  > Some other shapes are currently considered for special boundary shapes (sphere, cylinder...) and specifique loadings.
	
* **Boundary Conditions**: any rigid element can be used to apply boundary conditions. It is possible to impose velocity, force, or moment component by component. Some predefined systems with servo-control are also available for complex loading conditions (e.g., loading cycles or controlled pressure).

  > [!NOTE]
  > The possibility of applying tri-periodic loading to an assembly is implemented and currently in the testing phase.

* **Parallel Computation**: currently, an OpenMP optimization using compilation flags has been implemented. However, the computational speedup is relatively low. Typically, 8 cores are needed to halve the simulation time (for a dense system with a large number of elements).

* **Documentation**: there is little documentation, although efforts are being made to address this. For now, it is possible to generate the sphinxdoc documentation in your local folder.

<div style="text-align:center;">
<img src="./gif-doc/gen_doc-rockable.gif" width="90%"/>
</div>

## How to install

<div style="text-align:center;">
<img src="./gif-doc/install-script-rockable.gif" width="90%"/>
</div>

Using your OS package manager (yum, apt, brew etc) you will maybe need to install several package before compiling: `glfw3`, `opengl`,`freeglut`, `libpng2` (optionnal).

If you are lucky, the compilation is as simple as:

```sh
bash install_rockable.sh
```

The compilation is done in the `BUILD` directory and the binaries go to the `INSTALL` directory

Then, if needed, you can manage compilation options (profiling with MATools, full fetch of dependencies, see compilation, prepro compilation etc) using ccmake:

```sh
cd BUILD
ccmake .
# Set up options then `c`, then `e`, then `g`
cmake ..
make -j
make install 
```

## How to run 

Before runing rockable you will need to source rockable environnement to add the INSTALL directory to your standard binaries PATH:

```sh
source Env_rockable.sh
``` 

Then you can launch Rockable everywhere using:

```sh
rockable my_input_file.txt
```

