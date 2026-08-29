.. _drivingSystem:

Driving the boundaries
======================

In Rockable there is no special "wall" object: a boundary is an ordinary body
that happens to be driven. The bodies that can be driven are the **first**
``nDriven`` entries of the ``Particles`` list, and what is done to them is
described in a file named ``drivingSystem.txt``.

.. code-block:: text
   :caption: input.txt

   nDriven 6
   Particles 2006
   ...

.. contents::
   :local:
   :depth: 2


The file ``drivingSystem.txt``
------------------------------

The file sits in the folder where the computation runs. A token starting with
``/``, ``#`` or ``!`` discards the rest of the line.

.. note::

   If the file does not exist, the first ``nDriven`` bodies simply do not move.
   This is the usual way of building a container of fixed walls: no file is
   needed at all.

.. important::

   Unlike the rest of the setup, ``drivingSystem.txt`` is **re-read every time
   a conf-file is dumped**, that is every ``interConf``. A loading can therefore
   be changed while the computation runs, without restarting it: edit the file,
   and the change is taken into account at the next dump. This is also why the
   driving is never written into the conf-files.


Imposing a component: ``Control``
---------------------------------

The elementary entry drives **one component** of one body:

.. code-block:: text

   Control <type> <bodyNumber> <value>

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - ``<type>``
     - ``<value>``
     - Meaning
   * - ``_x_Vel_``, ``_y_Vel_``, ``_z_Vel_``
     - 1 real
     - Imposed velocity component.
   * - ``_xrot_Vel_``, ``_yrot_Vel_``, ``_zrot_Vel_``
     - 1 real
     - Imposed angular velocity component.
   * - ``_x_For_``, ``_y_For_``, ``_z_For_``
     - 1 real
     - Imposed force component.
   * - ``_xrot_Mom_``, ``_yrot_Mom_``, ``_zrot_Mom_``
     - 1 real
     - Imposed moment component.
   * - ``_xyzrot_Vel_``
     - 3 reals
     - Imposed angular velocity **vector**.
   * - ``_xyzrot_Mom_``
     - 3 reals
     - Imposed moment **vector**.

.. code-block:: text
   :caption: drivingSystem.txt

   # a floor moving up slowly
   Control _y_Vel_ 0 0.01

   # a piston pushed with a constant force
   Control _y_For_ 1 -250.0

   # a drum rotating about the z axis
   Control _xyzrot_Vel_ 2 0 0 1.5

A component that is **not** controlled stays free: the body is integrated along
that degree of freedom like any other one. A body listed in ``drivingSystem.txt``
with a single ``_y_Vel_`` control therefore falls under gravity along :math:`x`
and :math:`z` while its vertical motion is imposed. To keep a wall completely
fixed, either give it no control at all, or control all of its components.

.. note::

   Velocity-driven and force-driven components are not integrated in the same
   way. A force-driven component is integrated like a free body, with the
   imposed force added to the resultant; a velocity-driven component is simply
   translated at the imposed velocity, which makes it exactly kinematic. See
   :ref:`IntegrationSchemes`.

.. tip::

   ``TrackBody`` and ``TrackRockfall`` subtract the imposed force from the
   recorded resultant, so that the column holds the force actually transmitted
   by the contacts. This is how a wall is used as a force sensor.


Automatic controls: ``Servo``
-----------------------------

A ``Servo`` recomputes the values of its controls at every time step, from the
current state of the sample. It is the way to impose a **stress** rather than a
force, since the force to apply depends on the current area of the wall.

.. code-block:: text

   Servo <name> <parameters>

.. warning::

   Only **one** servo can be active: the servo function is a single callback,
   so the last ``Servo`` entry read replaces any previous one. Several
   ``Control`` entries, on the other hand, coexist without any problem.

The available servos are described in :ref:`Servo-controllers`.


Driving a periodic cell
-----------------------

When Rockable is compiled with ``ROCKABLE_ENABLE_PERIODIC``, the same file also
drives the periodic cell, through the ``PeriodicLoading`` keyword. There are no
walls in that case, and the loading is applied to the cell matrix itself. See
:ref:`periodicBoundaryConditions`.


A worked example: a box of six walls
------------------------------------

The most common setup is a cuboidal box whose six walls are the first six
bodies. With ``generator``, the walls are produced by
``generatePacking:wallBox``, which emits them in the order

.. code-block:: text

   0: Xmin   1: Xmax   2: Ymin   3: Ymax   4: Zmin   5: Zmax

Fixed box, particles simply poured in:

.. code-block:: text
   :caption: drivingSystem.txt (nothing to write)

   # no file needed: the 6 walls do not move

Oedometric compression, the top wall coming down at a constant velocity:

.. code-block:: text
   :caption: drivingSystem.txt

   Control _y_Vel_ 3 -0.01

Same, but pushing with a constant force instead:

.. code-block:: text
   :caption: drivingSystem.txt

   Control _y_For_ 3 -1000.0

Isotropic compression at a controlled **pressure**, which requires a servo
because the force must follow the changing wall areas:

.. code-block:: text
   :caption: drivingSystem.txt

   Servo tritriIsostaticCompression 0 1 2 3 4 5
   1000
