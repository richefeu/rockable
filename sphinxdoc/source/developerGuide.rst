.. _developerGuide:

Extending Rockable
==================

Rockable is meant to be extended by reading and modifying it, not by scripting
it from the outside. The design is deliberate: there is no Lua or Python layer,
so a new interaction model is written in C++, compiled in, and becomes a keyword
of the conf-file like any other.

Five families of components are pluggable through the same mechanism: force
laws, body forces, data extractors, pre-processing commands, and
post-processors.

.. contents::
   :local:
   :depth: 2


The factory mechanism
---------------------

Each family has an abstract base class, and the concrete classes register
themselves in a factory keyed by their **class name**. The registration happens
in ``Rockable::ExplicitRegistrations``:

.. code-block:: c++

   REGISTRER_BASE_DERIVED(ForceLaw, Avalanche);

The macro stringifies the class name, so the keyword written in the conf-file is
exactly the name of the C++ class. This is why ``forceLaw Avalanches`` does not
work while ``forceLaw Avalanche`` does: the factory is looked up by that exact
string.

.. code-block:: c++

   ForceLaw* FL = Factory<ForceLaw>::Instance()->Create(lawName);

When the lookup fails, Rockable warns and falls back to a default rather than
stopping, so a misspelt name produces a warning and a silently different
simulation. It is worth reading the first lines of the log.


The keyword parser
------------------

The conf-file is parsed by a ``kwParser``, a map from a keyword to a lambda that
reads the rest of the entry from the stream:

.. code-block:: c++

   parser.kwMap["tmax"] = __GET__(conf, tmax);

   parser.kwMap["density"] = __DO__(conf) {
     size_t grp;
     double density;
     conf >> grp >> density;
     properties.set(idDensity, grp, density);
   };

Anything the map does not know produces an ``unknown token`` message and is
skipped. Two consequences are worth keeping in mind:

.. warning::

   The map is keyed by string, so **two components registering the same keyword
   silently overwrite each other**, the last one registered winning. When adding
   a component, check that its keyword is not already taken.

.. note::

   Parsing is token-based, not line-based. The arguments of a keyword may be
   spread over several lines, which is what makes the long ``Servo`` and
   packing entries readable.


Adding a force law
------------------

1. Create ``src/ForceLaws/ForceLaw_MyLaw.hpp`` and ``.cpp``, deriving from
   ``ForceLaw``.

2. Implement ``init()``, which declares the parameters the law reads:

   .. code-block:: c++

      void MyLaw::init() {
        box->idKnContact = box->dataTable.add("knContact");
        box->idMuContact = box->dataTable.add("muContact");
      }

   ``dataTable.add`` both allocates the group-pair table and makes the keyword
   readable from the conf-file. A parameter that is not declared here cannot be
   set.

3. Implement ``computeInteraction(Interaction& I)``. It returns ``true`` when
   the interaction is active, that is when it carries a non-zero force, and
   ``false`` when it should be ignored:

   .. code-block:: c++

      bool MyLaw::computeInteraction(Interaction& I) {
        if (I.dn > 0.0) {          // no overlap, no contact
          I.fn = 0.0;
          I.ft.reset();
          I.mom.reset();
          return false;
        }
        ...
        return true;
      }

4. Register the class in ``Rockable::ExplicitRegistrations``, and add the file to
   ``CMakeLists.txt``.

.. important::

   ``computeInteraction`` is called from a parallel loop. It may write to the
   interaction it is given, but writing to a body, or to any shared structure,
   needs the protection that ``parallel_mode`` provides. Look at how the
   existing laws accumulate forces before adding your own.

.. tip::

   Start from ``ForceLaw_Default.cpp`` and change what you need.
   ``ForceLaw_geoVisc.cpp`` is a good example of a minimal departure: it is the
   default law with a viscous friction instead of an elastic one.


Adding a data extractor
-----------------------

Derive from ``DataExtractor`` and implement four methods:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Method
     - Role
   * - ``read(istream&)``
     - Read the parameters from ``dataExtractors.txt``, open the output file,
       and fill ``docString`` and ``columnDoc``.
   * - ``init()``
     - Anything that needs the configuration to be loaded and the driving
       system to be read.
   * - ``exec()``
     - Computations to perform with a period of ``nstep``.
   * - ``record()``
     - Write one line, with a period of ``nrec``.

.. code-block:: c++

   void MyExtractor::read(std::istream& is) {
     is >> ibody >> filename >> nrec;
     if (box->isInteractive() == false) recordFile.open(filename.c_str());
     nstep = std::numeric_limits<int>::max();

     docString << "ibody = " << ibody;
     columnDoc.clear();
     columnDoc.push_back("Time");
     columnDoc.push_back("what the second column holds");
   }

.. note::

   The ``isInteractive`` test matters: it prevents the extractor from opening,
   and truncating, its output file when the conf-file is merely being looked at
   with ``see`` or ``conftovtk``.

   Filling ``columnDoc`` is not optional in practice: it is what
   ``extractedDataDoc.txt`` is built from, and an output file whose columns are
   undocumented is of little use six months later.

Set ``nstep`` to a huge value when the extractor has nothing to compute between
two records, as all the current ones do.


Adding a pre-processing command
-------------------------------

Derive from ``PreproCommand`` and implement two methods:

.. code-block:: c++

   void myCommand::addCommand() {
     box->parser.kwMap["myCommand"] = [this](std::istream& conf) {
       conf >> this->someParameter;
       exec();
     };
   }

   void myCommand::exec() {
     // act on box->Particles, box->Interfaces, ...
   }

Then add the name to the ``commands`` array at the end of
``Rockable::initParser``, which is what instantiates the command and calls
``addCommand`` on it.

.. warning::

   The keyword written in ``addCommand`` must match the class name used in the
   array. This is precisely where a copy/paste error is easy to make and hard to
   notice: a wrong keyword there does not fail, it overwrites another command's
   entry and makes both behave unexpectedly.

.. note::

   A pre-processing command runs while the conf-file is being parsed, so it sees
   only what has been read **before** it. It must therefore be placed after the
   data it acts on, which is why they are documented as belonging at the end of
   the file.


Adding a post-processor
-----------------------

Derive from ``PostProcessor``, and implement ``read``, ``init``, ``exec`` and
``end``. The parameters are read with a local ``kwParser``, which lets a
post-processor have its own small keyword set:

.. code-block:: c++

   void MyPostProcessor::read(std::istream& is) {
     kwParser parser;
     parser.kwMap["Volume"] = __GET__(istr, Volume);
     parser.parse(is);
   }

``exec()`` is called once per processed conf-file, and ``end()`` after the last
one, which is where a summary or a distribution is written.


Where things live
-----------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Directory
     - Contents
   * - ``src/Core``
     - ``Rockable``, ``Particle``, ``Shape``, ``Interaction``,
       ``DrivingSystem``, ``PeriodicCell``: the heart of the code.
   * - ``src/Apps``
     - ``run.cpp`` (the ``rockable`` binary), ``see.cpp``, ``postpro.cpp``,
       ``conftovtk.cpp``, ``DearMySeer``, ``DearMyShape``.
   * - ``src/ForceLaws``
     - One file per interaction model.
   * - ``src/BodyForces``
     - One file per body force.
   * - ``src/DataExtractors``
     - One file per on-the-fly measurement.
   * - ``src/PreproCommands``
     - One file per pre-processing command.
   * - ``src/PostProcessors``
     - One file per post-processor.
   * - ``src/ProcessingTools``
     - Reusable analyses: clusters, broken sub-clusters, interaction groups,
       mass range, solid fraction probe.
   * - ``deps``
     - Dependencies fetched by CMake, ``toofus`` in particular.

The library ``toofus`` provides most of the utilities used throughout:
``vec3r``, ``quat``, ``mat9r``, ``AABB``, ``OBB``, ``kwParser``, ``DataTable``,
``Factory``, ``Tempo``, and the logging.


Testing a change
----------------

Build the regression tests with ``-DROCKABLE_USE_TESTING=ON``. The comparison
itself is available from the command line:

.. code-block:: sh

   rockable input.txt -n conf10 -r reference/conf10

which runs the simulation and then compares the two conf-files line by line from
their ``Particles`` section onwards, exiting with an error if they differ.

.. tip::

   A conf-file whose first token is ``redirection`` holds a path and a file
   name, and loads that file instead. It is how the tests share one sample
   between several cases without duplicating it.
