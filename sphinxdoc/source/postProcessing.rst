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


``SampleState``
---------------

The state of the whole sample, one line per processed configuration, in the
file ``SampleState.txt``: the time, the sizes and the volume of the sample, the
solid fraction, the coordination number, the number of contacts, and the stress
tensor computed from the contact forces.

.. code-block:: text
   :caption: commands.txt

   firstConf 0
   lastConf 100
   stepConf 1

   PostProcessor SampleState
   WallBox 0 1 2 3 4 5

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``WallBox <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>``
     - The numbers of the six walls of a box. The sizes of the sample are the
       distances between their inner faces (the position of each wall, minus
       or plus its Minkowski radius times its homothety).
   * - ``Volume <double>``
     - A volume imposed for every configuration (the sizes are then written
       as 0, or taken from the box or the cell if there is one).

With a periodic cell, nothing has to be given: the volume is
:math:`\det \mathbf{h}`, and the sizes are the lengths of the three vectors of
the cell.

The columns of ``SampleState.txt`` are:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Content
   * - 1
     - time
   * - 2-4
     - sizes :math:`L_x`, :math:`L_y`, :math:`L_z` of the sample
   * - 5
     - volume :math:`V`
   * - 6
     - solid fraction: the volume of the free particles (the volume of their
       shape times the cube of their homothety) divided by :math:`V`
   * - 7
     - coordination number :math:`Z`: twice the number of contacts between
       free bodies plus the number of contacts with driven bodies, divided by
       the number of free bodies. Two bodies touching at several points
       (vertex-face, edge-edge...) make **one** contact.
   * - 8
     - number of contacts (pairs of bodies)
   * - 9, 10
     - mean stress :math:`p = \operatorname{tr}\boldsymbol{\sigma}/3` and von
       Mises stress :math:`q`
   * - 11-16
     - :math:`\sigma_{xx}`, :math:`\sigma_{yy}`, :math:`\sigma_{zz}`,
       :math:`\sigma_{xy}`, :math:`\sigma_{xz}`, :math:`\sigma_{yz}`

The stress is

.. math::

   \boldsymbol{\sigma} = \frac{1}{V} \sum_{p} \sum_{c} \underline{F}_c \otimes
   \left( \underline{x}_p - \underline{x}_c \right)

where the first sum runs over the free bodies, the second over their contacts,
:math:`\underline{F}_c` is the force received by the body :math:`p` at the
contact point :math:`\underline{x}_c`. It is **positive in compression**, and
only its symmetric part is written. It does not include the kinetic
contribution of the velocity fluctuations, negligible in a quasi-static state.
With a periodic cell, it is the opposite of the stress that Rockable computes
during the computation to drive the cell (which is positive in tension).

.. note::

   With walls, the contact forces with the walls are included through the
   free bodies that touch them: the stress is that of the sample. In a box
   driven by a servo-control, :math:`\sigma_{xx}`, :math:`\sigma_{yy}` and
   :math:`\sigma_{zz}` match the imposed pressure once the sample is at
   equilibrium.


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
