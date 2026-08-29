.. _dataExtractors:

Data extractors
===============

Data extractors are probes evaluated **during** a computation. They write plain
text files that can be plotted directly, and they avoid the usual alternative of
saving hundreds of conf-files only to post-process a single scalar afterwards.

.. contents::
   :local:
   :depth: 2


The file ``dataExtractors.txt``
-------------------------------

The extractors are declared in a file named ``dataExtractors.txt``, placed in
the folder where the computation runs. The file is read once, at start-up, right
after the driving system. Each entry starts with the name of an extractor,
followed by its own parameters:

.. code-block:: text
   :caption: dataExtractors.txt

   # one line per extractor
   MeanVelocity     meanVelocity.txt 100
   TrackBody     12 body12.txt       100
   dnStat           dn.txt           500

A token starting with ``/``, ``#`` or ``!`` discards the rest of the line.

.. note::

   The extractors are **not** saved into the conf-files. The file is re-read
   from the running folder every time a computation is started, so restarting
   from ``conf27`` keeps the same probes, and editing the file between two runs
   is enough to change them.

.. warning::

   The extractors are ignored when the conf-file is opened by an interactive
   tool (``see``, ``seer``, ``conftovtk``, ``postpro``). Only ``rockable``
   evaluates them.


Common parameters
-----------------

Every extractor ends its parameter list with the same two entries:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Parameter
     - Description
   * - ``<filename>``
     - Name of the text file to write. It is opened at start-up and flushed
       after each record, so it can be plotted while the computation runs.
   * - ``<nrec>``
     - Number of **time steps** between two records. This is a step count, not
       a duration: with ``dt 1e-6``, ``nrec 1000`` records every millisecond of
       simulated time.

The first column of every output file is the time :math:`t`.

.. tip::

   A file named ``extractedDataDoc.txt`` is generated at start-up whenever at
   least one extractor is declared. It documents, extractor by extractor, what
   each column of each output file contains. It is the authoritative reference
   when a column ordering is not obvious.


Available extractors
--------------------

``MeanVelocity`` (*string*) **filename** (*int*) **nrec**
    Mean of the velocity magnitude :math:`\Vert \underline{v} \Vert` over
    **all** the bodies, driven ones included.

    Columns: ``t``, mean velocity.

``TrackBody`` (*int*) **ibody** (*string*) **filename** (*int*) **nrec**
    Full kinematics and loading of a single body.

    Columns: ``t``, then ``pos`` (3), ``vel`` (3), ``Q`` (4, as ``w x y z``),
    ``vrot`` (3), the resultant force (3) and the resultant moment (3).

    .. note::

       When a component of the body is driven in *force*, the imposed force is
       subtracted from the recorded resultant, so that the column holds the
       force actually carried by the contacts.

``TrackRockfall`` (*int*) **ibody** (*double*) **vStop** (*double*) **wStop** (*string*) **filename** (*int*) **nrec**
    Same columns as ``TrackBody``, but it also **stops the computation** when
    the tracked block has come to rest, which is what a rockfall study needs:
    there is no point integrating a block that has stopped.

    The run is stopped as soon as
    :math:`\Vert \underline{v} \Vert < v_\text{stop}` **and**
    :math:`\Vert \underline{\omega} \Vert < w_\text{stop}`, or as soon as the
    block leaves the bounding box of the driven bodies.

    .. note::

       The stop test is inhibited before :math:`t_\text{min} = v_\text{stop}/g`
       with :math:`g = 9.81`. This grace period is computed automatically, and
       prevents a block released at rest from stopping the computation at the
       very first step, before gravity has had time to accelerate it.

``ClusterAABB`` (*int*) **icluster** (*string*) **filename** (*int*) **nrec**
    Axis-aligned bounding box of one cluster, useful to follow the swelling or
    the spreading of a fragmenting body.

    Columns: ``t``, ``aabb.min`` (3), ``aabb.max`` (3).

``dnStat`` (*string*) **filename** (*int*) **nrec**
    Statistics of the normal distances :math:`d_n` over all the interactions.
    It is the quickest way to check that a simulation is not over-compressed
    and that the time step is small enough.

    Columns: ``t``, ``dnMin``, ``dnMax``, ``dnMean``, the mean over the negative
    values only, the mean over the positive values only, the count of negative
    values, the count of positive values, and the position (3) where
    ``dnMin`` occurs.

``DuoBalance`` (*int*) **i** (*int*) **j** (*string*) **filename** (*int*) **nrec**
    Composition of the interaction between two given bodies, by sub-interaction
    type.

    Columns: ``t``, number of vertex-vertex, vertex-edge, vertex-face and
    edge-edge sub-interactions, then their sum weighted by the contact
    partnership.

``TrackDamage`` (*string*) **filename** (*int*) **nrec**
    Damage of the sample, defined as the broken interface area divided by the
    interface area present initially.

    Columns: ``t``, damage in :math:`[0, 1]`.

    .. note::

       Only meaningful with a force law that breaks interfaces
       (``StickedLinks``, ``BCM``) and once the interfaces have been created by
       a sticking pre-processing command.


Declaring an extractor inside a conf-file
-----------------------------------------

A ``DataExtractor`` keyword is also accepted inside a conf-file:

.. code-block:: text

   DataExtractor MeanVelocity meanVelocity.txt 100

.. warning::

   This form is kept for backward compatibility only, and it emits a warning.
   It has a real drawback: the line is **not** written back into the
   conf-files that the computation saves, so restarting from a dump silently
   loses the extractor. Use ``dataExtractors.txt`` instead.


Plotting the results
--------------------

The output files are column-formatted text, directly usable by ``gnuplot``:

.. code-block:: text

   plot 'meanVelocity.txt' u 1:2 w l t 'mean velocity'
   plot 'body12.txt' u 1:3 w l t 'height of body 12'

The same holds for the three files that ``rockable`` always writes
(``perf.txt``, ``kineticEnergy.txt`` and ``staticBalance.txt``); the generated
script ``checkplots.txt`` plots them in one go:

.. code-block:: sh

   gnuplot checkplots.txt
