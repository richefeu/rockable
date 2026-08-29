.. _shapeSurvey:

Application ``shapeSurvey``
===========================

``shapeSurvey`` is an OpenGL browser for a shape library. It answers the two
questions that come up whenever a shape file has been written by hand or
produced by a converter: *does this shape look like what I meant?* and *are its
mass properties right?*

.. code-block:: sh

   shapeSurvey shapes.txt

It is also the tool that **pre-computes** a library. A shape whose
``preCompDone`` is ``n`` has its volume, inertia and bounding box recomputed
every time a simulation starts; running the computation once here and saving the
result removes that cost, and lets you check the values before they are used.

.. contents::
   :local:
   :depth: 2


Browsing the library
--------------------

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Key
     - Action
   * - ``+`` / ``-``
     - Next / previous shape. A shape whose ``preCompDone`` is ``n`` gets its
       bounding box fitted on the fly when it is displayed.
   * - ``h``
     - Show the help.
   * - ``e``
     - Print the extents of the bounding box in the terminal.
   * - ``a`` / ``A``
     - Decrease / increase the transparency.
   * - ``b``
     - Background colours on/off.
   * - ``w`` / ``W``
     - Roll the camera about the viewing axis.
   * - ``q``
     - Quit.

The mouse rotates, pans and zooms as in ``see``.


Computing the mass properties
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Key
     - Action
   * - ``c``
     - Compute the mass properties of the current shape, and set its
       ``preCompDone`` to ``y``.
   * - ``C``
     - The same for **every** shape of the library that is still marked ``n``.
   * - ``*``
     - Reset ``preCompDone`` of the current shape to ``n``, to force a
       recomputation.
   * - ``N`` / ``n``
     - Multiply / divide by ten the number of Monte-Carlo samples
       (``MCnstep``), between :math:`10^3` and :math:`10^8`.
   * - ``d``
     - Clean every shape, removing the duplicated entities.

The volume and the inertia of a sphero-polyhedron have no closed form in the
general case, so they are estimated by Monte-Carlo sampling. The accuracy
therefore depends on ``MCnstep``, and the estimated relative error on the volume
is printed with the result.

.. tip::

   Increase ``MCnstep`` with ``N`` until the printed error is small enough for
   your purpose, then press ``C`` and ``s``. For a shape used as a driven wall
   the default is plenty; for a grain whose mass drives the dynamics, it is
   worth a few more samples.

.. note::

   A shape marked ``preCompDone y`` is **trusted**: Rockable uses the volume and
   the inertia written in the file without checking them. This is efficient, and
   it is also how a wrong value gets used silently, so it is worth looking at
   the numbers once.


Fitting the bounding box
------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Key
     - Action
   * - ``o``
     - Cycle through the fitting strategies: covariance, minimum volume,
       axis-aligned, imposed axis. The chosen one is stored as
       ``fibObbOption``.
   * - ``t``
     - Build the OBB-tree of the current shape.
   * - ``k`` / ``K``
     - Show one level less / one level more of the tree.

The oriented bounding box is what the neighbour detection uses first, so a
poorly fitted box costs time on every step. The four strategies differ in what
they optimise: the covariance fit is fast and usually good, the minimum-volume
fit is tighter on elongated shapes, and the axis-aligned one is only relevant
for shapes already aligned with the axes.

.. warning::

   The OBB-tree, and the ``OBBtreeLevel`` keyword of the shape files, are
   deprecated. The tree can still be built and displayed here, but the
   contact detection is driven by ``AddOrRemoveInteractions``, whose
   ``OBBtree`` option builds what it needs on its own.


Saving
------

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - Key
     - Action
   * - ``s``
     - Save the library back to **the file it was read from**.
   * - ``p``
     - Export a sample: one particle per shape, laid out so that the library can
       be opened directly in ``see``.

.. warning::

   ``s`` overwrites the input file in place, without asking. Keep a copy of a
   library you care about before pressing it.


Typical workflow after a conversion
-----------------------------------

.. code-block:: sh

   stl2shape -i mesh.stl -r 0.01 -c    # -> mesh.shp, preCompDone n
   shapeSurvey mesh.shp

Then, in the window: browse with ``+`` to check the geometry, press ``C`` to
compute the mass properties of every shape, and ``s`` to save. The library is
then ready for a simulation, and no shape will be recomputed at start-up.
