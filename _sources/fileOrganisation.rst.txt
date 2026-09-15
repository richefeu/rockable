.. _File-organisation:

Organisation of the files
=========================


📁 ``examples``
---------------

Directory containing example simulations or processings. The example name **helloworld** can be a starting point but the reste of the folder is a kind of sand-box the developers use to test their features.


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

Directory containing non-regression tests for the project. It would be nice to add more tests.
