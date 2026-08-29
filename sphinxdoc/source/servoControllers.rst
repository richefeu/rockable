.. _Servo-controllers:

Available servo-controllers
===========================

A servo-controller is a function evaluated at every time step that rewrites the
values of the controls it owns. It is declared in ``drivingSystem.txt`` with the
``Servo`` keyword (see :ref:`drivingSystem`), and it is what makes
**stress-controlled** and **time-dependent** loadings possible: a stress cannot
be imposed directly, because the force to apply depends on the current area of
the wall, which changes as the sample deforms.

.. warning::

   Only one servo can be active at a time. The servo function is stored as a
   single callback, so a second ``Servo`` entry silently replaces the first.

.. contents::
   :local:
   :depth: 2


The ``tritri`` family
---------------------

These servos drive a cuboidal cell made of **six walls**, and they all begin
with the same six body numbers:

.. code-block:: text

   Servo tritri<Something> <idXmin> <idXmax> <idYmin> <idYmax> <idZmin> <idZmax>
   <parameters...>

They install twelve controls at once, in a fixed pattern: the walls at the
**minimum** positions are velocity-driven, and the walls at the **maximum**
positions are force-driven.

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Wall
     - Control installed
     - Role
   * - ``idXmin``, ``idYmin``, ``idZmin``
     - ``_x_Vel_``, ``_y_Vel_``, ``_z_Vel_``
     - Held fixed (velocity imposed at zero), they define the origin.
   * - ``idXmax``, ``idYmax``, ``idZmax``
     - ``_x_For_``, ``_y_For_``, ``_z_For_``
     - Driven by a force recomputed at every step from the target stress.

At each step, the servo measures the current inner dimensions of the cell from
the wall positions, corrected by their Minkowski radii:

.. math::

   d_X = \left| X_\text{max} - R_\text{max} - \left( X_\text{min} + R_\text{min} \right) \right|

and likewise for :math:`d_Y` and :math:`d_Z`. The areas of the faces follow,
:math:`S_X = d_Y d_Z`, :math:`S_Y = d_X d_Z`, :math:`S_Z = d_X d_Y`, and the
force applied to a wall is the target stress times the corresponding area.

.. note::

   Because the Minkowski radii are subtracted, the dimensions used are those of
   the **free volume** available to the grains, not the distance between the
   wall centres.


``tritriIsostaticCompression`` (*double*) **pressure**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Isotropic compression. The same pressure is applied on the three ``max`` walls,
the three ``min`` walls staying fixed:

.. math::

   f_X = -p\, S_X \qquad f_Y = -p\, S_Y \qquad f_Z = -p\, S_Z

.. code-block:: text
   :caption: drivingSystem.txt

   Servo tritriIsostaticCompression 0 1 2 3 4 5
   1000

This is the standard way of preparing a dense sample: compress isotropically
until the kinetic energy has dropped, then use the resulting conf-file as the
initial state of the real test.


``tritriBiaxialCompression`` (*double*) **pressure** (*double*) **velocity**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compression along :math:`y` at an imposed velocity, with a confining pressure
kept constant on the lateral directions :math:`x` and :math:`z`.

The control of the ``Ymax`` wall, installed as a force by the common pattern, is
converted back into a velocity control set to :math:`-v`, while the servo keeps
recomputing the lateral forces:

.. math::

   f_X = -p\, S_X \qquad f_Z = -p\, S_Z

.. code-block:: text
   :caption: drivingSystem.txt

   Servo tritriBiaxialCompression 0 1 2 3 4 5
   1000 0.01


``tritriCustom`` (*int*) **type** (*double*) **value** ... (twelve values)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The general form: each of the six walls is given its own type and value, in the
order ``Xmin``, ``Xmax``, ``Ymin``, ``Ymax``, ``Zmin``, ``Zmax``.

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - ``type``
     - Meaning of ``value``
   * - ``0``
     - Imposed velocity of the wall.
   * - ``1``
     - Imposed **stress**, converted to a force with the current wall area at
       every step.

.. code-block:: text
   :caption: drivingSystem.txt

   Servo tritriCustom 0 1 2 3 4 5
   0 0.0
   1 1000.0
   0 0.0
   0 -0.01
   0 0.0
   1 1000.0

This example holds the three ``min`` walls fixed, confines :math:`x` and
:math:`z` at 1000 Pa, and moves the top wall down at 0.01 m/s, which reproduces
``tritriBiaxialCompression`` and shows how to depart from it.

.. note::

   The values of the ``max`` walls are applied with a reversed sign, so that a
   positive stress always means compression whatever the wall.


``tritriLodeAngle`` (*double*) **pressure** (*double*) **LodeAngle** (*double*) **sigRate**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A stress path at constant Lode angle, used to explore the deviatoric plane
rather than a single triaxial path. Starting from an isotropic state at
``pressure``, the stress increment grows linearly in time,
:math:`\Delta\sigma = \dot{\sigma} t`, and is distributed over the three
directions as

.. math::

   \sigma_X = p + \Delta\sigma \qquad
   \sigma_Y = p + a \Delta\sigma \qquad
   \sigma_Z = p + b \Delta\sigma

with

.. math::

   a = \frac{3 \tan\theta_L - \sqrt{3}}{2\sqrt{3}} \qquad b = -(1 + a)

where :math:`\theta_L` is the Lode angle, in degrees, restricted to
:math:`[0°, 60°]`.

.. warning::

   This servo has not been fully tested. It also assumes that the computation
   starts at :math:`t = 0`, since :math:`\Delta\sigma` is computed from the
   absolute time: restarting from a conf-file at :math:`t \neq 0` would resume
   with an already non-zero stress increment.


Shakers
-------

Shakers impose an oscillating motion on a single body, along a given direction.
They install three velocity controls (``_x_Vel_``, ``_y_Vel_``, ``_z_Vel_``) on
that body and rewrite them at every step. The direction is normalised
automatically, so it does not need to be a unit vector.

.. note::

   What is imposed is the **velocity**, taken as the exact derivative of the
   intended motion. The body therefore oscillates about the position it had
   when the shaking started, with amplitude :math:`A`, rather than about a
   position given in the file.


``shaker`` (*int*) **body** (*vec3r*) **dir** (*double*) **A** (*double*) **freq**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sinusoidal oscillation of amplitude :math:`A` and frequency :math:`f`:

.. math::

   \underline{v}(t) = A\, \omega \cos(\omega t)\, \underline{d}
   \qquad \omega = 2 \pi f

.. code-block:: text
   :caption: drivingSystem.txt

   Servo shaker 0  0 1 0  0.001 50


``triangle_shaker`` (*int*) **body** (*vec3r*) **dir** (*double*) **A** (*double*) **freq**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same parameters, but the velocity is a **square** wave, so the displacement is a
triangular wave: the body moves at the constant speed :math:`4Af` and reverses
twice per period. The acceleration is zero except at the reversals, which avoids
the continuously varying inertial forcing of a sine.

``sawtooth_shaker`` (*int*) **body** (*vec3r*) **dir** (*double*) **A** (*double*) **freq** (*double*) **t_ini**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``triangle_shaker`` with a settable phase origin ``t_ini``, so that the
oscillation can be made to start at a chosen time rather than at :math:`t = 0`.
This is what makes it usable after a deposition stage.


Ramps
-----

``ramp`` (*string*) **type** (*int*) **body** (*double*) **valueBegin** (*double*) **valueEnd** (*double*) **tBegin** (*double*) **tEnd**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Drives a single control, of any of the types listed in :ref:`drivingSystem`,
along a linear ramp:

.. math::

   v(t) = \begin{cases}
     v_\text{begin} & t \leq t_\text{begin} \\
     v_\text{begin} + \dfrac{v_\text{end} - v_\text{begin}}{t_\text{end} - t_\text{begin}} \left( t - t_\text{begin} \right) & t_\text{begin} < t < t_\text{end} \\
     v_\text{end} & t \geq t_\text{end}
   \end{cases}

.. code-block:: text
   :caption: drivingSystem.txt

   # push harder and harder between t = 0.1 s and t = 0.6 s
   Servo ramp _y_For_ 3 0.0 -5000.0 0.1 0.6

This is the way to avoid the shock of a load applied abruptly at the first step.

.. tip::

   A ramp on a *material parameter* rather than on a boundary is a different
   mechanism: the ``Tempo`` keyword of the conf-file, described in
   :ref:`syntaxConf`, which can ramp a friction coefficient or the numerical
   damping over time.
