
Format of shape-files
=====================

Definition of a shape
---------------------

A shape-file is a library of shapes, which is in fact a file with a number of shapes defined between brackets (starting with ``<`` and ending with ``>``). To run a simulation, we need that all required shapes are defined in a single shape-file. The keywords are given here after:

- ``name`` (*string*) **shapeName**
  The name of the shape
  
- ``radius`` (*double*) **MinskowskiRadius**
  The radius of the rounded edges (so-called Minskowski radius) of the shape
  
- ``preCompDone`` ``y`` or ``n``
  If yes, the computation of mass properties will not be run at the beginning of a computation 
  in Rockable or in the program shapeSurvey
  
- ``OBBtreeLevel`` (*int*) **number**
  set the number of levels for building an OBB-tree
  
- ``nv`` (*int*) **numberOfVertices**

  | Then repeat **numberOfVertices** times: 
  | (*double*) **x** (*double*) **y** (*double*) **z**

- ``ne`` (*int*) **numberOfEdges**

  | Then repeat **numberOfEdges** times: 
  | (*int*) **from** (*int*) **to**

- ``nf`` (*int*) **numberOfFaces**

  | Then repeat **numberOfFaces** times: 
  | (*int*) **numberOfVertices** [(*int*) **V1** (*int*) **V2** ...]
  
- ``obb.extent`` (*double*) **extent1** (*double*) **extent2** (*double*) **extent3**
  The 3 components of the extents in the 3 directions defined by ``obb.e1``, ``obb.e2`` and ``obb.e3``
  
- ``obb.e1`` (*vec3r*) **e1**
  Extent unit vector in the first direction
  
- ``obb.e2`` (*vec3r*) **e2**
  Extent unit vector in the second direction
  
- ``obb.e3`` (*vec3r*) **e3**
  Extent unit vector in the third direction
  
- ``obb.center`` (*vec3r*) **center**
  Position of the OBB center relative to the mass center of the shape, given in the framework (principal direction of   
  inertia) of the shape
  
- ``volume`` (*double*) **V**
  Volume of the shape
  
- ``I/m`` (*double*) **I1/m** (*double*) **I2/m** (*double*) **I3/m**
  Eigenvalues of the moment of inertia divided by the mass. It supposes that the positions of the vertices 
  are given in the eigen-frame of the shape

- ``position`` (*vec3r*) **pos**
  A position that can be stored for processing purpose

- ``orientation`` (*quat*) **Q**
  As the position defined previouly, an orientation can also be stored for processing purpose 

As an example, this is how a cube with unit side can be defined (radius of 0.1):

.. code-block:: text

   <
   name Cube_r0.1
   radius 0.1
   preCompDone y

   nv 8
   0.4 0.4 -0.4
   -0.4 0.4 -0.4
   -0.4 -0.4 -0.4
   0.4 -0.4 -0.4
   0.4 0.4 0.4
   -0.4 0.4 0.4
   -0.4 -0.4 0.4
   0.4 -0.4 0.4

   ne 12
   0 1
   1 2
   2 3
   3 0
   4 5
   5 6
   6 7
   7 4
   0 4
   1 5
   2 6
   3 7
   
   nf 6
   4 0 1 2 3
   4 4 5 6 7
   4 0 1 5 4
   4 2 3 7 6
   4 1 2 6 5
   4 0 4 7 3
   
   obb.extent 0.5 0.5 0.5
   obb.e1 1 0 0
   obb.e2 0 1 0
   obb.e3 0 0 1
   obb.center 0 0 0

   volume 0.975587
   I/m 0.166667 0.166667 0.166667
   >




