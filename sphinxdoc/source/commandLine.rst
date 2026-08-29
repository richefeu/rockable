.. _commandLine:

Command line and build options
==============================

.. contents::
   :local:
   :depth: 2


``rockable``
------------

.. code-block:: sh

   rockable [conf-file] [options]

The positional argument is the conf-file to load; it defaults to ``conf0``.
Everything happens in the **current folder**: the companion files are looked up
there, and the dumps are written there.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``-j``, ``--nbThreads`` (*int*)
     - Number of OpenMP threads, 1 by default. The parallel *strategy* itself is
       a conf-file keyword, ``parallel_mode``.
   * - ``-v``, ``--verbose`` (*int*)
     - Verbosity, from 0 to 6: ``off``, ``critical``, ``err``, ``warn``,
       ``info`` (the default), ``debug``, ``trace``.
   * - ``-c``, ``--clean``
     - Delete ``conf*``, ``kineticEnergy.txt``, ``perf.txt``,
       ``staticBalance.txt`` and ``checkplots.txt`` from the current folder,
       then exit.
   * - ``-b``, ``--banner``
     - Print the banner, which lists the version, the git tag and the
       compilation options actually used, then exit.
   * - ``-n`` (*string*) ``-r`` (*string*)
     - Non-regression check. After the run, compare the two conf-files line by
       line from their ``Particles`` section onwards, and exit with an error if
       they differ.
   * - ``-h``, ``--help``
     - The usual help.
   * - ``--version``
     - The git tag the binary was built from.

.. code-block:: sh

   rockable input.txt -j 8 -v 5

.. warning::

   ``rockable -c`` deletes files without asking for a confirmation. Run it in
   the folder of the simulation you really mean to clean.

.. tip::

   Any dump is a valid input file, so ``rockable conf27 -j 8`` resumes the
   computation from configuration 27. The companion files are re-read from the
   folder, which is how a loading can be changed mid-course.


The other executables
---------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Command
     - Description
   * - ``see [conf-file]``
     - The OpenGL viewer. ``-t <file>`` loads a trajectory file, ``-v <int>``
       sets the verbosity of the core. See :ref:`Visualisation`.
   * - ``seer``
     - *Dear My Seer*, the Dear ImGui viewer that supersedes ``see``.
   * - ``conftovtk``
     - Convert **every** ``conf*`` of the folder into VTK files for ParaView.
   * - ``postpro <commandFile>``
     - Run a post-processor over a range of dumps. See :ref:`postProcessing`.
   * - ``generator <commandFile>``
     - Build shape files, input files and packings. See the ``generator`` page.
   * - ``shapeSurvey <shapeFile>``
     - Inspect and pre-compute a shape library. See :ref:`shapeSurvey`.
   * - ``sweepable [script]``
     - Generate a tree of input files for a parameter sweep.
   * - ``stl2shape -i <f.stl> -r <R>``
     - Convert a binary STL file into a shape. See :ref:`converters`.


Building
--------

.. code-block:: sh

   git clone https://github.com/richefeu/rockable.git
   cd rockable
   sh install_rockable.sh

The build happens in ``BUILD`` and the binaries are installed into ``INSTALL``.
Depending on the system, a few packages may be needed beforehand: ``glfw3``,
``opengl``, ``freeglut``, and optionally ``libpng`` for PNG screenshots.

To change the options:

.. code-block:: sh

   cd BUILD
   ccmake .
   # set the options, then c, then e, then g
   cmake ..
   make -j
   make install

To make the binaries reachable from the shell:

.. code-block:: sh

   source add_install_to_path.sh

.. important::

   The script must be **sourced**, not executed. A script run in a subshell
   cannot change the ``PATH`` of the shell that launched it.


Compilation options
-------------------

Options that change the physics, or the set of keywords understood. All default
to ``OFF``.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Option
     - Effect
   * - ``ROCKABLE_ENABLE_PERIODIC``
     - Tri-periodic cell: the ``usePeriodicCell``, ``h``, ``vh``, ``ah``,
       ``mh``, ``dh`` keywords and the ``PeriodicLoading`` drivings. See
       :ref:`periodicBoundaryConditions`.
   * - ``ROCKABLE_ENABLE_SOFT_PARTICLES``
     - Homogeneous straining of the particles, keyword ``useSoftParticles``.
       See :ref:`softParticles`.
   * - ``ROCKABLE_ENABLE_BOUNDARY``
     - The special boundary shapes, ``Ball`` and ``Cylinder``.
   * - ``ROCKABLE_USE_FT_CORR``
     - Objectivity correction of the tangential forces under large rotations.
   * - ``ROCKABLE_ENABLE_PROFILING``
     - Time profiling of the main routines, reported into ``perf.txt``.

.. warning::

   ``ROCKABLE_ENABLE_BOUNDARY`` compiles the ``Ball`` and ``Cylinder``
   boundaries, but their ``read`` methods are still empty stubs in the source.
   Enabling the option therefore does not yet give a usable feature.

Options that select what gets compiled:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Option
     - Default
     - Builds
   * - ``ROCKABLE_COMPILE_SEE``
     - ``ON``
     - ``see``
   * - ``ROCKABLE_COMPILE_SEER``
     - ``ON``
     - ``seer`` (Dear My Seer)
   * - ``ROCKABLE_COMPILE_CONF2VTK``
     - ``ON``
     - ``conftovtk``
   * - ``ROCKABLE_COMPILE_POSTPRO``
     - ``ON``
     - ``postpro``
   * - ``ROCKABLE_COMPILE_PREPRO``
     - ``ON``
     - the pre-processing tools
   * - ``ROCKABLE_USE_TESTING``
     - ``OFF``
     - the regression tests

.. tip::

   ``rockable -b`` prints which of these were actually enabled. It is the
   quickest way to explain a keyword that appears to do nothing.


Files of a simulation folder
----------------------------

.. code-block:: text

   mySimulation/
     input.txt              # the conf-file you write
     shapes.txt             # shape library, named by shapeFile
     shapes.rmsh            # optional skin meshes, used by see
     drivingSystem.txt      # optional driving and servo
     dataExtractors.txt     # optional measurements
     see.json               # optional view settings
     probe.txt              # optional measurement box
     ------------ produced by the run ------------
     conf0 conf1 conf2 ...  # the dumps
     extractedDataDoc.txt   # what each extractor column holds
     perf.txt               # time, efficiency, profiling
     kineticEnergy.txt      # time, translational and rotational energy
     staticBalance.txt      # time, equilibrium indicators
     checkplots.txt         # gnuplot script for the three files above

.. tip::

   ``gnuplot checkplots.txt`` plots the performance, the kinetic energy and the
   static balance in one go. Watching the kinetic energy fall is the usual way
   of deciding that a sample has reached equilibrium.
