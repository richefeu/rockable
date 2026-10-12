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

Computes, for each processed configuration, the tensorial moment of each
particle,

.. math::

   \mathbf{M}_i = \sum_c \underline{f}_c \otimes \underline{x}_c

where the sum runs over the interactions of the particle, :math:`\underline{f}_c`
is the **opposite** of the force it receives (so that a compression is
positive) and
:math:`\underline{x}_c` the position of the contact point. For a particle in
equilibrium, this is the usual sum of the forces times the branch vectors. The
driven bodies are skipped.

The results are written into ``ParticleStress-conf<iconf>.txt``, one line per
particle: the 9 components of :math:`\mathbf{M}_i` (they are **not** divided
by the volume of the particle), then four values that separate an axial part,
along the axis :math:`z` of the body (:math:`\underline{u}`), from a radial
part: :math:`a = \underline{u} \cdot \mathbf{M}_i \underline{u}`,
:math:`r = \sqrt{r_1^2 + r_2^2}` with :math:`r_k = \underline{t}_k \cdot \mathbf{M}_i \underline{t}_k`
for two unit vectors :math:`\underline{t}_1, \underline{t}_2` perpendicular to
:math:`\underline{u}`, then the trace of the axial part and the trace of the
radial part. The stress of the whole sample, :math:`\sum_i \mathbf{M}_i / V`, is
printed on the terminal.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``Volume <double>``
     - Total volume used to normalise the stress of the sample. This is the
       default value, used for every configuration that is not listed below.
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

.. note::

   For the stress of the soft particles themselves, which Rockable computes
   during the run, use the data extractors ``TrackSoftParticle`` or
   ``TrackStress`` (see :ref:`dataExtractors`).


``ClusterGranulo``
------------------

Distribution of the sizes of the sub-clusters, for fragmentation studies: as
the interfaces break, an initial cluster splits into sub-clusters, and this
post-processor counts them.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``SievingSizes <nb>``
     - Followed by ``nb`` values, one per class.

.. code-block:: text
   :caption: commands.txt

   firstConf 0
   lastConf 50
   stepConf 10

   PostProcessor ClusterGranulo
   SievingSizes 5
   1 2 3 4 5

The results are written into ``ClusterGranulo-conf<iconf>.txt``, one line per
class: the value given in ``SievingSizes``, and the fraction (in number) of the
sub-clusters in that class.

.. warning::

   Despite their name, the values of ``SievingSizes`` are not sizes that the
   sub-clusters are compared with: they are only written as labels. The class
   of a sub-cluster is its **number of particles** minus one (a single particle
   is in the first class, two particles in the second, ...), and the larger
   sub-clusters all go into the last class. It is therefore a distribution of
   the number of particles per sub-cluster, not a grain-size distribution.

.. note::

   The sub-clusters are identified from the interfaces that are still intact,
   so this analysis is only meaningful together with a force law that breaks
   interfaces (``StickedLinks`` or ``BCM``). The file of a configuration
   without any sub-cluster is left empty.


Post-processing without ``postpro``
-----------------------------------

The conf-files are plain text, and their format is documented in
:ref:`syntaxConf`. For a one-off analysis, reading them with a short script is
often quicker than writing a new post-processor. The ``Particles`` and
``Interactions`` blocks are column-formatted, and the number of entries is given
on the line that opens each block.

Writing a new post-processor in C++ is described in :ref:`developerGuide`.
