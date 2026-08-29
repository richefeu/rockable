.. _postProcessing:

Post-processing (``postpro``)
=============================

Where the :ref:`data extractors <dataExtractors>` measure things *during* a
computation, ``postpro`` works *afterwards*, by re-reading the conf-files that
were saved. This is the right tool for quantities that are too expensive to
evaluate at every step, or that were not anticipated when the run was launched.

.. code-block:: sh

   postpro commands.txt

.. note::

   ``postpro`` is built when the CMake option ``ROCKABLE_COMPILE_POSTPRO`` is
   on, which is the default.


The command file
----------------

The command file selects a range of conf-files and one post-processor:

.. code-block:: text
   :caption: commands.txt

   firstConf 0
   lastConf 100
   stepConf 5

   PostProcessor ParticleStress
   Volume 1.0

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Keyword
     - Description
   * - ``firstConf <int>``
     - Number of the first conf-file to process.
   * - ``lastConf <int>``
     - Number of the last conf-file to process (inclusive).
   * - ``stepConf <int>``
     - Increment between two processed conf-files.
   * - ``PostProcessor <name>``
     - Selects the post-processor. The keywords that follow, until the end of
       the file, are its own parameters.

The conf-files are loaded one after the other, from ``./conf<first>`` to
``./conf<last>``, and the post-processor is executed on each of them.

.. warning::

   ``postpro`` loads the configurations in **interactive mode**. The data
   extractors declared in ``dataExtractors.txt`` are therefore not run again,
   and nothing is written back into the conf-files.

.. important::

   Only **one** post-processor can be used per command file. Run ``postpro``
   twice, with two command files, to obtain two different analyses.


``ParticleStress``
------------------

Computes the stress tensor of the sample from the contact forces and the branch
vectors, for each processed configuration.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``Volume <double>``
     - Total volume used to normalise the stress. This is the default value,
       used for every configuration that is not listed below.
   * - ``ConfVolumes <nb>``
     - Followed by ``nb`` lines ``<iconf> <volume>``, giving the volume of
       specific configurations.

The volume has to be provided because Rockable has no general way of knowing
the volume actually occupied by a sample: it depends on the boundaries, which
may be walls, a periodic cell, or nothing at all. When the sample is compacted
during the run, its volume changes from one configuration to the next, and
``ConfVolumes`` is the way to follow it:

.. code-block:: text
   :caption: commands.txt

   firstConf 0
   lastConf 2
   stepConf 1

   PostProcessor ParticleStress
   Volume 1.0
   ConfVolumes 3
   0 0.0012
   1 0.00121
   2 0.00123

.. tip::

   The six bounds of a box of driven walls are printed by ``see`` with the key
   ``x``, which gives a quick way of building this list for a triaxial cell.


``ClusterGranulo``
------------------

Grain-size distribution of the clusters, obtained by sieving. It is meant for
fragmentation studies: as the interfaces break, an initial cluster splits into
sub-clusters, and this post-processor counts how many pass through each sieve.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``SievingSizes <nb>``
     - Followed by ``nb`` sieve sizes, given as plain numbers on the following
       lines.

.. code-block:: text
   :caption: commands.txt

   firstConf 0
   lastConf 50
   stepConf 10

   PostProcessor ClusterGranulo
   SievingSizes 5
   0.001 0.002 0.005 0.01 0.02

.. note::

   The sub-clusters are identified from the interfaces that are still intact,
   so this analysis is only meaningful together with a force law that breaks
   interfaces (``StickedLinks`` or ``BCM``).


Post-processing without ``postpro``
-----------------------------------

The conf-files are plain text, and their format is documented in
:ref:`syntaxConf`. For a one-off analysis, reading them with a short script is
often quicker than writing a new post-processor. The ``Particles`` and
``Interactions`` blocks are column-formatted, and the number of entries is given
on the line that opens each block.

Writing a new post-processor in C++ is described in :ref:`developerGuide`.
