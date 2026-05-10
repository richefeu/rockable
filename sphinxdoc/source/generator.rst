Application `generator`
=======================

``generator`` reads commands from a text file and produces Rockable
configuration files (shape files, input files, packings). It is the standard
pre-processing tool for setting up simulations.

.. code-block:: sh

   generator <command-file>

.. note::

   Leave at least one blank line at the end of the command file to avoid
   parsing issues.

.. contents::
   :local:
   :depth: 2


Example
-------

The following script generates a shape file and a full Rockable input file
with a random close packing of cuboids:

.. code-block:: text

   # --- shapes ---
   open shapes.txt
   generateShape:cuboid Cuboid 0.0005 0.03 0.01 0.01
   generateShape:xyz_walls 0.25 0.25 0.25 0.0005
   close

   # --- input file ---
   open input.txt
   print Rockable 21-08-2022
   print #### time #############################
   print t 0
   print tmax 2.0
   print dt 5e-6
   print interVerlet 0.001
   print interConf 0.01
   print
   print #### options ##########################
   print forceLaw Default
   print Integrator Beeman
   print
   print #### proximity of the particles #######
   print DVerlet 0.005
   print dVerlet 0.0025
   print
   print #### properties for particles #########
   print density 0 2700
   print density 1 2700
   print
   print #### interaction parameters ###########
   print knContact 0 1 1e6
   print en2Contact 0 1 0.02
   print ktContact 0 1 1e6
   print muContact 0 1 0.6
   print krContact 0 1 0.0
   print murContact 0 1 0.0
   print
   print #### the system #######################
   print iconf 0
   print shapeFile shapes.txt
   print nDriven 0
   Particles:placeHolder

   generatePacking:RandomClosePacking
   Cuboid
   CUBOID Y
   0.2 0.2 0.3
   0.03 0.01 0.01
   1.0 1.0
   2000
   0.9
   0
   0

   close


Command Reference
-----------------

File output
^^^^^^^^^^^

``open <filename>``
   Open ``filename`` for output. All subsequent ``print``, ``compute``, and
   packing commands write to this file until ``close`` is called. If no file
   is open, output goes to the console.

``close``
   Close the current output file. If ``Particles:placeHolder`` was used, the
   final particle count is written back into the file at that position before
   closing.

Printing and arithmetic
^^^^^^^^^^^^^^^^^^^^^^^

``print <text>``
   Write ``text`` to the current output, trimming leading and trailing spaces,
   and append a newline. Useful for writing raw lines into shape or
   configuration files.

``print> <text>``
   Same as ``print`` but without the trailing newline. The next output
   statement continues on the same line.

   Example — writing a computed contact stiffness on one line:

   .. code-block:: text

      print> knContact 0 0
      <compute 500000 * 2

   produces: ``knContact 0 0 1000000``

``compute <label> <expression>``
   Evaluate a mathematical expression and write ``label result`` to the output.

   Example: ``compute Particles 5^3 + 6`` produces ``Particles 131``

``<compute <expression>``
   Evaluate an expression and write only the result (preceded by a space),
   without a label. Intended to be combined with ``print>`` on the previous
   line (see example above).

Particle counting
^^^^^^^^^^^^^^^^^

``Particles:placeHolder``
   Reserve a slot in the output file where the total particle count will be
   written when ``close`` is called. The count is updated automatically by
   every packing command.

``incrementNoParticles <n>``
   Manually increase the tracked particle count by the integer ``n``. Use this
   when adding particles through ``print`` rather than through a packing
   command.

Transformations
^^^^^^^^^^^^^^^

A global transformation is maintained throughout script execution and applied
by packing commands when placing particles.

``transformation:add_translation <x> <y> <z>``
   Add a translation to the current transformation.

``transformation:add_rotation <x> <y> <z> <angle_deg>``
   Add a rotation around the axis ``(x, y, z)`` by ``angle_deg`` degrees.

``transformation:reset``
   Reset the transformation to identity (clears both translation and rotation).

``transformation:reset_translation``
   Reset only the translation component.

``transformation:reset_rotation``
   Reset only the rotation component.


Shape Generation
----------------

``generateShape:sphere <name> <radius>``
   Generate a spherical shape.

``generateShape:cube <name> <radius> <sideSize>``
   Generate a cubic shape of side length ``sideSize`` with corner radius
   ``radius``.

``generateShape:cuboid <name> <radius> <sx> <sy> <sz>``
   Generate a cuboid shape with side lengths ``(sx, sy, sz)`` and corner
   radius ``radius``.

``generateShape:pyramid3 <name> <radius> <sideSize>``
   Generate a triangular pyramid shape.

``generateShape:thin_cylinder <name> <Rin> <Rout> <H> <nbSectors>``
   Generate a thin cylindrical shell with inner radius ``Rin``, outer radius
   ``Rout``, height ``H``, and ``nbSectors`` angular sectors.

``generateShape:cuboid_container <name> <radius> <sx> <sy> <sz> <hasTop> <hasBottom>``
   Generate an open or closed cuboid container. ``hasTop`` and ``hasBottom``
   are integers (``0`` = absent, ``1`` = present).

``generateShape:rectangle_xz <name> <radius> <side_x> <side_z>``
   Generate a rectangular plate in the *xz* plane.

``generateShape:rectangle_xy <name> <radius> <side_x> <side_y>``
   Generate a rectangular plate in the *xy* plane.

``generateShape:rhombicuboctahedron <name> <radius> <sx> <sy> <sz>``
   Generate a rhombicuboctahedron shape.

``generateShape:xyz_walls <sx> <sy> <sz> <Rw>``
   Generate three flat wall shapes (one per axis) sized ``(sx, sy, sz)`` with
   corner radius ``Rw``.


Packing Generation
------------------

``addParticle <name> <group> <cluster> <homothety> <px> <py> <pz> <qs> <qx> <qy> <qz>``
   Place a single particle at position ``(px, py, pz)`` with quaternion
   orientation ``(qs, qx, qy, qz)``, homothety scale ``homothety``, assigned
   to ``group`` and ``cluster``.

``generatePacking:wallBox <group> <LX> <LY> <LZ> <Rw>``
   Generate six flat wall particles forming a closed box of dimensions
   ``LX × LY × LZ`` with wall corner radius ``Rw``, all assigned to
   ``group``.

``generatePacking:grid <name> <ox> <oy> <oz> <bx> <by> <bz> <nx> <ny> <nz> <group> <cluster> <homothety> <randQ>``
   Fill a box of size ``(bx, by, bz)`` starting at origin ``(ox, oy, oz)``
   with a regular ``nx × ny × nz`` grid of particles of shape ``name``.
   Set ``randQ`` to ``1`` for random orientations.

``generatePacking:grid_clust <name> <ox> <oy> <oz> <bx> <by> <bz> <nx> <ny> <nz> <group> <cluster> <homothety> <randQ>``
   Same as ``generatePacking:grid`` but assigns incrementing cluster indices
   to each particle.

``generatePacking:RandomClosePacking <name> <boxShape> <direction> <bx> <by> <bz> <xOBB> <yOBB> <zOBB> <hmin> <hmax> <nbTarget> <solidFraction> <group> <cluster>``
   Generate a random close packing inside a box of size ``(bx, by, bz)``.

   .. list-table::
      :header-rows: 1
      :widths: 25 75

      * - Parameter
        - Description
      * - ``boxShape``
        - Container geometry: ``CUBOID`` or ``CYLINDER``
      * - ``direction``
        - Packing growth direction: ``X``, ``Y``, or ``Z``
      * - ``xOBB`` ``yOBB`` ``zOBB``
        - Bounding box dimensions of the particle shape
      * - ``hmin`` ``hmax``
        - Homothety range (size scatter)
      * - ``nbTarget``
        - Target number of particles
      * - ``solidFraction``
        - Target solid fraction (0–1)
      * - ``group`` ``cluster``
        - Group and cluster indices assigned to all particles