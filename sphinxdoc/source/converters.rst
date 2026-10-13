.. _converters:

Converters
==========

Four converters bridge Rockable's shape format and the usual geometry formats.
They live in ``prepro/converters``, and are built by CMake with
``ROCKABLE_COMPILE_PREPRO``.

.. contents::
   :local:
   :depth: 2


``stl2shape``: from a triangulated mesh
---------------------------------------

Converts a **binary** STL mesh into a single shape, by taking the mesh as the
skeleton of a sphero-polyhedron and dilating it by a Minkowski radius.

.. code-block:: sh

   stl2shape -r 0.01 -c file.stl

The shape is written in the file ``mesh.shp`` of the current folder (an existing
``mesh.shp`` is overwritten), and is named after the STL file, without its
folder (``file.stl``). The options ``-o`` and ``-n`` change both.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Option
     - Description
   * - *file.stl* (last argument)
     - Input STL file. Required, and it must be the **binary** flavour of STL.
   * - ``-r``, ``--radius`` (*double*)
     - Minkowski radius of the resulting shape. Required.
   * - ``-s``, ``--scaleFactor`` (*double*)
     - Rescale the object by a given factor.
   * - ``-m``, ``--maxLength`` (*double*)
     - Rescale the object so that its largest dimension, its sieving size, takes
       the given value.
   * - ``-z``, ``--scaleRadius``
     - Rescale the radius along with the object, instead of keeping it as given.
   * - ``-c``, ``--clean``
     - Remove the duplicated edges.
   * - ``-n``, ``--name`` (*string*)
     - Name of the shape. By default, the name of the STL file, without its
       folder.
   * - ``-o``, ``--output`` (*string*)
     - Name of the shape file written, ``mesh.shp`` by default.

.. important::

   Use ``-c`` almost always. An STL file stores each triangle independently, so
   every edge appears twice and every vertex several times. Without cleaning,
   the shape carries duplicated sub-elements, which costs time at every contact
   detection and can produce forces counted twice.

.. warning::

   The cost of a shape grows with its number of vertices, edges and faces, and
   so does the cost of the contact detection. A mesh exported at full resolution
   from a CAD tool or a scanner easily has thousands of triangles, which is far
   more than a discrete element simulation needs. Decimate the mesh **before**
   converting, and consider ``AddOrRemoveInteractions OBBtree`` for the complex
   shapes that remain.

The converted shape has ``preCompDone n``: run ``shapeSurvey`` on it to compute
its mass properties and save them (see :ref:`shapeSurvey`).


``tess2shape``: from a Neper tessellation
-----------------------------------------

Converts a `Neper <https://neper.info>`_ tessellation into a shape library
**and** the matching particle list. This is the usual route to a sample of
space-filling grains, of the kind used for the bonded models.

.. code-block:: sh

   tess2shape command.txt

The cells are shrunk by the Minkowski radius (each face moves inwards by the
radius), so that two neighbouring cells touch without overlapping, and can be
glued by ``stickBCM``. The mass properties of each cell are computed by a
Monte-Carlo integration of ``MCnstep`` points: most of the time of the
conversion, about 30 s for 1000 cells with ``20000``.

.. note::

   ``tif2rockable``, in the same folder, is still **not built by CMake** (its
   line is commented out in ``prepro/converters/CMakeLists.txt``).

The command file is a list of keywords:

.. code-block:: text
   :caption: command.txt

   tessFileName        tessel.tess
   inputFileName       input.txt
   shapesFileName      shapes.txt
   MinkowskiRadius     0.01
   ParticlesGroup      0
   ParticlesCluster    0
   ParticlesHomothety  1.0
   MCnstep             50000

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``tessFileName``
     - The ``.tess`` file produced by Neper.
   * - ``inputFileName``
     - Conf-file to write, holding the particle list.
   * - ``shapesFileName``
     - Shape library to write, one shape per cell of the tessellation.
   * - ``MinkowskiRadius``
     - Radius given to every generated shape.
   * - ``ParticlesGroup``, ``ParticlesCluster``, ``ParticlesHomothety``
     - Group, cluster and homothety given to every generated particle.
   * - ``MCnstep``
     - Number of Monte-Carlo samples used for the mass properties.

.. warning::

   The keyword is ``MinkowskiRadius``. Older command files may spell it
   ``MinskowskiRadius``, which is **not** recognised: the line is ignored and
   the radius silently keeps its default value.

Generating the tessellation itself is a Neper matter:

.. code-block:: sh

   neper -T -n 100 -domain "cube(1,1,1)" -o tessel

Neper places its domain at the origin; ``-transform "translate(x,y,z)"`` moves
it, for instance to surround a shape that will be used as a cookie cutter
(see ``examples/bunnyFun/voronoiBunny``).

.. warning::

   Neper 5.0.1 crashes with a domain read from a file,
   ``-domain "file(shape.obj)"``, even for a cube: a non-convex domain cannot
   be tessellated directly. The tutorial ``voronoiBunny`` tessellates a box,
   and keeps the cells whose centre is inside the shape.

.. tip::

   Neper is built from its sources (``https://github.com/neperfepx/neper``),
   with ``gsl`` and ``nlopt``. With GCC 16, which compiles C in C23 by default,
   the library Scotch that Neper includes needs
   ``-DCMAKE_C_FLAGS="-std=gnu17 -Wno-incompatible-pointer-types"``.

.. tip::

   Giving every particle the same ``cluster`` makes all the bonds *inner*, which
   is what a single fragmenting solid needs. Distinct clusters make them
   *outer*, for an assembly of separate grains. See :ref:`prePro`.


``stl2blocks``: a mesh cut into blocks
--------------------------------------

Cuts a closed STL mesh (binary or ascii) into blocks along families of parallel
planes, and writes the blocks as a shape library **and** the matching particle
list, ready to be glued by ``stickBCM``.

.. code-block:: sh

   stl2blocks command.txt

The cuts form a tree. The whole mesh (the block ``root``) is cut along one
direction into slabs, numbered ``0``, ``1``, ... along that direction; each slab
can then be cut along its own direction, with its own spacing of the planes,
and so on. A block that is not cut any more is a final block (or several, when
the planes leave several separate parts of the mesh in it).

The Minkowski radius ``R`` is given first, and the blocks are built so that,
once rounded by ``R``, they fill the mesh and touch each other:

- the faces of the mesh are moved inwards by ``R`` (the *core* of the mesh);
- each cutting plane removes a slice of thickness ``2R``: the faces of the two
  blocks on either side of it are ``2R`` apart.

.. code-block:: text
   :caption: command.txt

   stlFile          ../bunny.stl
   homothety        2
   translation      0 0.0515 0
   MinkowskiRadius  0.001

   cut root      1 0 0      spacing 0.02 0
   cut root.*    0 1 0.3    spacing 0.02 0.005
   cut root.*.*  0.2 -0.3 1 number 5
   cut root.3.*  1 1 1      random 2 3

   shapesFile       blocks.shp
   inputFile        blocks.txt
   ParticlesGroup   1
   ParticlesCluster 1
   MCnstep          10000

A line ``cut <path> <direction> <planes>`` cuts the blocks whose path matches
``<path>``: ``root``, ``root.2`` (the slab 2 of ``root``), ``root.2.0``... A
``*`` matches any slab number; when several lines match a block, the one with
the fewest ``*`` is used (the last one, if they have as many). Above,
``root.3.*`` is cut along ``1 1 1``, and the other slabs of the second level
along ``0.2 -0.3 1``. The direction need not be normalised. The planes are placed
in one of five ways, along the direction (a position ``p`` is the plane
``x.u = p``, ``u`` being the normalised direction):

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Planes
     - Description
   * - ``number <n>``
     - ``n`` slabs of the same thickness, over the extent of the block.
   * - ``spacing <s> <offset>``
     - The planes ``offset + k s``: the same planes for all the blocks cut with
       this line (they line up from one block to the next).
   * - ``positions <k> <p1> ... <pk>``
     - ``k`` planes at the given positions.
   * - ``fractions <k> <f1> ... <fk>``
     - ``k`` planes at the given fractions (between 0 and 1) of the extent of
       the block.
   * - ``random <k> <seed>``
     - ``k`` planes at random positions in the extent of the block.

Only the planes inside the block are used. The program prints the tree of cuts,
with the slabs kept (a slab with no part of the mesh is not kept, but the
numbers of the other slabs do not change), then one line per final block: its
name, its path, its number of vertices and faces, its volume.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``stlFile``
     - The STL file, binary or ascii.
   * - ``homothety``, ``translation``
     - Applied to the mesh before the cuts (and before ``MinkowskiRadius``,
       which is in the units of the result).
   * - ``MinkowskiRadius``
     - Radius of all the blocks (required).
   * - ``cut``
     - A cut (see above), one line per rule.
   * - ``shapesFile``, ``inputFile``
     - The shape library and the conf-file (the line ``shapeFile`` and the
       particles) to write (``blocks.shp`` and ``blocks.txt`` by default).
   * - ``shapeName``
     - The blocks are named ``<shapeName>_0``, ``<shapeName>_1``...
       (``block`` by default).
   * - ``ParticlesGroup``, ``ParticlesCluster``
     - Group and cluster of all the particles (``0`` by default). The same
       cluster for all makes all the bonds *inner* (see ``tess2shape``).
   * - ``MCnstep``
     - Number of Monte-Carlo points for the mass properties of each block
       (``20000`` by default).
   * - ``minVolume``
     - The parts whose core is smaller than this volume are removed (``0``).
   * - ``gluePoints``
     - ``1`` (default) or ``0``: the vertices added for ``stickBCM`` (below).

**The vertices for the glue.** ``stickBCM`` bonds the vertices of two blocks
that face each other, ``2R`` apart. Two blocks on either side of a plane are
cut by different later planes, and their vertices do not face each other: a
cubic grid of cuts gives four bonds per interface, an oblique tree of cuts
often one or two, and an interface with fewer than three bonds has no area
(it breaks at once). So the program adds to both blocks the corners of the
part of the plane they share: the vertices of each one that face the other
block, and the points where their sides cross. These vertices are in the faces
of the blocks, and do not change their shapes. The program also prints the area
of the faces of the cuts, to be compared with the area given by ``stickBCM``.

**Repair of the mesh.** The mesh is checked first: the triangles stuck on the
mesh by a single side (a *fin*) are removed, the triangles are oriented
outwards, and the holes are closed by fans of triangles. The Stanford bunny of
``examples/bunnyFun`` has one fin.

.. warning::

   The parts of the mesh thinner than ``2R`` cannot be represented. Moving the
   faces inwards by ``R`` is done node by node: on triangles smaller than
   ``R``, the nodes are moved less, so that no triangle is turned over (the
   program says how many); the core may then cross itself, and a face of a cut
   have a crossed outline (the program says so too). Choose ``R`` small
   compared with the details of the mesh, and with the spacing of the planes:
   a slab thinner than ``2R`` disappears, and a block whose slabs are all
   thinner than ``2R`` is not cut.

.. tip::

   Two parts of a block closer than ``2R`` to each other (a thin neck of the
   mesh, cut by the planes) would overlap once rounded: they are kept together,
   in one block. Small blocks need a small time step: the program gives the
   smallest one. For the bunny of ``examples/bunnyFun`` (homothety 2) cut as
   above, ``dt 2e-7`` is stable, ``1e-6`` is not. Even so, 29 of its 780
   interfaces (those with fewer than three bonds, where two blocks touch only
   along a sliver) break at the first step.

The program works on a copy of ``polyhTool.hpp`` of toofus (in its folder),
extended to polyhedra that are not convex: the cut by a plane gives faces with
several loops, holes, and concave outlines, triangulated by ear clipping.


``shape2mesh``: from a shape to a surface mesh
----------------------------------------------

The reverse direction: meshes the **exact** surface of a sphero-polyhedron, that
is the Minkowski sum of the skeleton with a ball. Rounded edges, corners,
concavities and open surfaces are all handled, and the vertex normals are the
analytic normals of the surface, so even a coarse mesh shades correctly.

.. code-block:: sh

   shape2mesh shapes.txt --rmsh

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Option
     - Description
   * - ``-f``, ``--format`` ``obj|ply``
     - Output format, ``obj`` by default.
   * - ``-e``, ``--epsilon`` (*double*)
     - Absolute sag tolerance, in the length unit of the shape file. Smaller
       means finer; the grid step follows :math:`h \sim \sqrt{8 \epsilon R}`.
       The default is 1% of the radius.
   * - ``-s``, ``--simplify`` [*degrees*]
     - Merge the coplanar triangles, with a default tolerance of 0.25 degree,
       and re-triangulate the flat regions from their boundary only.
   * - ``-r``, ``--rmsh``
     - Write a single ``<stem>.rmsh`` file holding every shape of the library.
   * - ``-o``, ``--outdir`` (*string*)
     - Output directory.

The mesh is expressed in the body frame of the shape, the frame stored in the
shape file.

The ``.rmsh`` companion
^^^^^^^^^^^^^^^^^^^^^^^

With ``-r``, all the shapes of a library are meshed into one ``.rmsh`` file.
When that file sits next to the shape file, ``see`` draws each particle with its
skin mesh instead of the overlapping spheres, tubes and slabs it normally uses.
The difference is purely visual, but on a rounded polyhedron it is the
difference between a figure that reads and one that does not.

.. note::

   The mesh is built in the body frame of the shape: its mass centre and its
   principal axes, where Rockable and ``see`` put a shape when they read it.
   A shape still marked ``preCompDone n``, such as the output of
   ``stl2shape``, is first pre-computed by ``shape2mesh`` (it says so, and it
   takes the time of a Monte-Carlo integration). The shape file itself is not
   modified: pre-compute it once with ``shapeSurvey shapes.shp -P shapes.shp``
   (see :ref:`shapeSurvey`) to save that time, here and at each start of a
   computation.

The skin mesh cannot follow the deformation of a soft particle (see
:ref:`softParticles`): the soft particles are always drawn with the overlapping
primitives, deformed.

.. code-block:: sh

   shape2mesh shapes.txt --rmsh       # -> shapes.rmsh
   shape2mesh shapes.txt --rmsh -s    # smaller, flat faces merged

Without an explicit ``-e``, the sag is chosen per shape, at 0.4% of its
bounding-box diagonal, so a large flat wall and a small detailed grain both get
a sensible mesh.

.. note::

   The simplification is safe on the flat parts: the offset faces of a
   sphero-polyhedron are exact planes, so merging their triangles does not
   change the surface, and the mesh stays watertight. The curved regions, the
   rounded edges and corners, are left untouched. On a rounded polyhedron the
   triangle count typically drops by about 60%; on a sphere, nothing changes.
