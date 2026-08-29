.. _converters:

Converters
==========

Three converters bridge Rockable's shape format and the usual geometry formats.
They live in ``prepro/converters`` and are built with
``ROCKABLE_COMPILE_PREPRO``.

.. contents::
   :local:
   :depth: 2


``stl2shape``: from a triangulated mesh
---------------------------------------

Converts a **binary** STL mesh into a single shape, by taking the mesh as the
skeleton of a sphero-polyhedron and dilating it by a Minkowski radius.

.. code-block:: sh

   stl2shape -i file.stl -r 0.01 -c

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Option
     - Description
   * - ``-i``, ``--input`` (*string*)
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

   The keyword is ``MinkowskiRadius``. Some example files shipped with the
   repository spell it ``MinskowskiRadius``, which is **not** recognised: the
   line is ignored and the radius silently keeps its default value.

Generating the tessellation itself is a Neper matter:

.. code-block:: sh

   neper -T -n 100 -domain "cube(1,1,1)" -o tessel

.. tip::

   Giving every particle the same ``cluster`` makes all the bonds *inner*, which
   is what a single fragmenting solid needs. Distinct clusters make them
   *outer*, for an assembly of separate grains. See :ref:`prePro`.


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
