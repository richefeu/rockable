.. _File-organisation:

Organisation of the files
=========================


📁 ``examples``
---------------

Directory containing example simulations or processings. The example named **helloworld** can be a starting point. The examples of **softParticles** run in a few seconds and each has a ``readme.md`` telling what to check (see :ref:`softParticles`). The rest of the folder is a kind of sandbox the developers use to test their features.


📁 ``prepro``
-------------

Some preprocessing applications.

- 📁 ``common``
  Shared between the editors: ``rockable.lang``, which describes the input-file
  language (keywords, types, documentation and snippets), and the header that
  reads it.

- 📁 ``confedit``
  Graphical tool (FLTK) for editing configurations.

- 📁 ``rockedit``
  Terminal editor for the same files, with the same colouring and documentation.
  No dependency, so it can be built and used over ssh.

- 📁 ``converters``
  Conversion applications.

  - 📁 ``shape2mesh``
    Builds the skin mesh (``.rmsh``) of the shapes, used by ``see`` to draw them.

  - 📁 ``stl2blocks``
    Cuts a closed STL mesh into blocks (shapes and particles), along planes.

  - 📁 ``stl2shape``
    Converter for STL files to shape format.

  - 📁 ``tess2shape``
    Converter for tessellation files to shape format.
    
  - 📁 ``tif2rockable``
    Converter for TIF files to a format usable for rock simulations.

- 📁 ``genesis``
  Directory for data generation tools.

  - 📁 ``cpptools``
    Tools written in C++ for data generation.

  - 📁 ``SpherePacker``
    Tool for sphere packing.

  - 📁 ``TubePacker``
    Tool for tube packing.

- 📁 ``shapesurvey``
  Tool for shape analysis.

- 📁 ``sweepable``
  Tool to run parametric studies (see :ref:`sweepable`).



📁 ``src``
----------

Main source directory containing the source code of ``Rockable``.

- 📁 ``Apps``
  Applications of the project: ``rockable`` for computation, ``see`` for visualisation, and other tools.

- 📁 ``BodyForces``
  Classes for body forces.

- 📁 ``Boundaries``
  Classes for boundary conditions.

- 📁 ``Core``
  Core classes of the project.

- 📁 ``DataExtractors``
  Classes for data extraction.

- 📁 ``ForceLaws``
  Classes for force laws.

- 📁 ``PostProcessors``
  Classes for post-processing.

- 📁 ``PreproCommands``
  Classes for preprocessing.

- 📁 ``ProcessingTools``
  Processing tools.

📁 ``test``
-----------

Directory containing the tests of the project: unit tests of the skin mesh, and
non-regression tests that run short simulations and compare their last
conf-file with a reference. They are built with ``-DROCKABLE_USE_TESTING=ON``
and run with ``ctest``; the file ``test/README.md`` explains how they work and
how to add one.


📁 ``sphinxdoc``
----------------

The sources of this documentation. The file ``sphinxdoc/README.md`` explains how
it is written and maintained.


📁 ``cheatsheets``
------------------

One-page reference sheets, written with Typst, for the conf-files, the companion
files, the physics, the command lines, ``see`` and the other tools.


📁 ``deps``
-----------

The libraries fetched by CMake at the configuration step (``toofus``, ``fmt``,
...). It is not part of the repository.
