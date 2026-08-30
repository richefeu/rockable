.. _periodicBoundaryConditions:

Periodic boundary conditions
============================

Periodic boundary conditions replace the walls by a repeating cell, which
removes the edge effects that a container inevitably introduces. This matters
whenever the quantity of interest is a bulk property rather than the response of
a finite sample, and it is what makes concurrent double-scale simulations
possible.

.. important::

   These keywords are only understood when Rockable has been compiled with

   .. code-block:: sh

      cmake .. -DROCKABLE_ENABLE_PERIODIC=ON

   Without that option the keywords are still accepted, but each of them only
   prints ``PERIODIC_NOT_ENABLED when Rockable was compiled`` and does nothing.
   ``rockable -b`` tells you which options a binary was built with.

.. warning::

   Periodic boundary conditions and the special boundaries
   (``ROCKABLE_ENABLE_BOUNDARY``) cannot be used together.

.. contents::
   :local:
   :depth: 2


The cell
--------

The cell is described by a :math:`3 \times 3` matrix :math:`\mathbf{h}` whose
**columns** are the three vectors spanning the cell. A cubic cell of side
:math:`L` is simply :math:`\mathbf{h} = L\,\mathbf{I}`, and a sheared cell has
off-diagonal terms.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Keyword
     - Description
   * - ``usePeriodicCell`` (*int*)
     - ``1`` activates the periodic cell, ``0`` disables it.
   * - ``h`` (*mat9r*)
     - The cell matrix, nine reals given row by row.
   * - ``vh`` (*mat9r*)
     - Velocity of the cell matrix, :math:`\dot{\mathbf{h}}`.
   * - ``ah`` (*mat9r*)
     - Acceleration of the cell matrix, :math:`\ddot{\mathbf{h}}`.
   * - ``mh`` (*double*)
     - Mass ratio used to build the inertia of the cell degrees of freedom.
   * - ``dh`` (*double*)
     - Numerical damping applied to the cell degrees of freedom.

The cell has its own dynamics: it is not a kinematic constraint but a set of
nine degrees of freedom integrated along with the particles, which is what
allows a **stress** to be imposed on it. ``mh`` sets how heavy those degrees of
freedom are, and therefore how fast the cell responds; ``dh`` damps their
oscillations.

.. tip::

   A cell that oscillates instead of converging under a constant pressure is
   the sign of an ``mh`` too small or a ``dh`` too small. Both are numerical
   parameters, to be calibrated like ``numericalDampingCoeff``.


Corrections and stress measurement
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Keyword
     - Description
   * - ``cellVelocityCorrection`` (*int*)
     - Removes the mean velocity drift of the sample.
   * - ``cellMomentumCorrection`` (*int*)
     - Removes the total momentum drift of the sample.
   * - ``useKineticStress`` (*int*)
     - Includes the kinetic (velocity fluctuation) term in the stress used to
       drive the cell.
   * - ``cellRelattice`` (*int*)
     - ``1`` keeps the cell on a short basis under a large shear. Off by
       default.
   * - ``cellDriveCauchy`` (*int*)
     - ``1`` makes a stress-driven component of the cell reach equilibrium on
       the corresponding component of the Cauchy stress. ``0``, the default, is
       the exact Parrinello-Rahman conjugate of :math:`\mathbf{h}`.

A fully periodic system has no boundary to anchor it, so nothing prevents the
whole sample from drifting: any residual momentum is conserved for ever. The two
correction flags remove that drift, and are normally left on.

``useKineticStress`` matters as soon as the sample is not quasi-static. The
stress of a granular assembly is the sum of a contact term and a kinetic term;
for a slow compression the second is negligible, for a rapid flow it is not.

``cellRelattice`` is only needed for large shear strains. A periodic system is
defined by its lattice, not by the box used to draw it: any
:math:`\mathbf{h}' = \mathbf{h}\,\mathbf{M}` with :math:`\mathbf{M}` an
integer matrix of determinant one describes the same system. Under a shear the
driven off-diagonal term of :math:`\mathbf{h}` grows without bound, the cell
leans further and further, and the perpendicular width of the cell shrinks as

.. math::

   w = \frac{h_{xx}}{\sqrt{1 + \gamma^2}}

Once :math:`w` falls below twice the interaction range, the nearest image found
by rounding the reduced coordinates is no longer the nearest one in space, and
contacts are silently missed. With ``cellRelattice 1``, whenever an off-diagonal
term exceeds half of the diagonal it leans on, the cell is re-expressed on a
shorter basis of the same lattice. Positions, velocities and accelerations
follow the change of basis, so nothing moves and the tangential history of the
contacts is preserved.

.. tip::

   For a cell of about seven grain diameters, the width criterion is only
   reached around :math:`\gamma \approx 1.4`; below that the option changes
   nothing. Switch it on when a run is meant to accumulate a shear strain of
   order one or more.

``cellDriveCauchy`` matters as soon as the cell is sheared. The force that
drives a stress-driven component of :math:`\mathbf{h}` is, by default and as in
the Parrinello-Rahman formulation, the exact conjugate of that component,

.. math::

   F_{ij} = V \sum_k h^{-1}_{ik} \, \Delta\Sigma_{kj}

The row of :math:`\mathbf{h}^{-1}` appearing here is a reciprocal cell vector,
which leans as soon as the cell does. For a simple shear of amplitude
:math:`\gamma`, row :math:`x` reads :math:`(1/h_{xx},\, -\gamma/h_{xx},\, 0)`
and the :math:`xx` component therefore settles on

.. math::

   \Sigma_{xx} - \gamma\, \Sigma_{xy} = -p
   \qquad\text{instead of}\qquad
   \Sigma_{xx} = -p

Stretching :math:`\mathbf{a}_1` along :math:`x` while :math:`\mathbf{a}_2` is
tilted is not a pure normal strain, so the conjugate force legitimately mixes the
two stress components. The trouble is that :math:`\gamma` is a property of the
basis, not of the lattice: the very same physical state written on another basis
of the same lattice is then held at another pressure. Measured on a shear at
:math:`p = 100` with the re-lattice on, :math:`\Sigma_{xx}` swept from
:math:`81.6` to :math:`120.0` in step with :math:`\gamma \in [-0.5, 0.5]`,
while :math:`\Sigma_{yy}` and :math:`\Sigma_{zz}`, whose rows of
:math:`\mathbf{h}^{-1}` do not lean, stayed at :math:`100.0`.

This is a property of the loading rather than an implementation error, and it
has a practical consequence: the mean pressure moves away from the setpoint as
the shear accumulates, so a friction coefficient has to be normalised by the
measured pressure and not by the setpoint. ``cellRelattice`` bounds the excursion
instead of letting it grow, and averaging over one full period of the resulting
sawtooth cancels it.

``cellDriveCauchy 1`` departs from that formulation. Only the diagonal term of
:math:`\mathbf{h}^{-1}` is kept, so a stress-driven component equilibrates on
its own component of the Cauchy stress whatever the shape of the cell, and two
bases of one lattice drive the cell identically. The driving force is then no
longer the variational conjugate of :math:`\mathbf{h}`; it is the one that
imposes what the loading says it imposes. Use it when the setpoint has to be
held literally; leave it off to stay on the reference formulation.

.. note::

   Kinematics are stored in the conf-files in **real** coordinates, and
   converted to the reduced coordinates of the cell internally when the file is
   read. A conf-file of a periodic simulation therefore remains readable, and
   ``see`` displays it correctly.


Driving the cell
----------------

The loading is set in ``drivingSystem.txt`` (see :ref:`drivingSystem`), with the
``PeriodicLoading`` keyword instead of the ``Control`` and ``Servo`` entries used
for walls.

``PeriodicLoading IsotropicCompression`` (*double*) **pressure**
    The three diagonal components of the cell are stress-driven at
    :math:`-p`, and every shear component is held at zero velocity. This is the
    periodic equivalent of ``tritriIsostaticCompression``, and the usual way of
    preparing a dense periodic sample.

``PeriodicLoading TriaxialCompression`` (*string*) **X|Y|Z** (*double*) **pressure** (*double*) **strainRate**
    The named direction is strain-driven at a constant rate, the two others
    being kept at ``pressure``. The imposed velocity is recomputed at every step
    as :math:`\dot{\varepsilon}\, h_{ii}`, so the **strain rate** stays
    constant as the cell shrinks, rather than the velocity.

``PeriodicLoading SimpleShearDeformable`` (*string*) **XY|XZ|YX|YZ|ZX|ZY** (*double*) **pressure** (*double*) **shearRate** [``FixedTransverse``]
    Shear at a constant rate on the named component, the three normal components
    being stress-driven at ``pressure``. The two letters name the component of
    :math:`\mathbf{h}` that is sheared; the normal direction used to convert the
    rate into a velocity follows from it.

    With the optional ``FixedTransverse``, only the dimension carrying the
    normal stress is stress-driven; the two others are held at their initial
    value, as in a shear box. The sample is then free to dilate only across the
    shear planes, which is what a laboratory direct-shear test does, and the two
    transverse dimensions can no longer respond to the stress at all.

.. code-block:: text
   :caption: drivingSystem.txt

   PeriodicLoading IsotropicCompression 100

.. note::

   ``PeriodicLoading`` replaces the walls entirely: ``nDriven`` should be 0, and
   there is no wall body in the ``Particles`` list.


Visualising a periodic sample
-----------------------------

``see`` draws the cell as a box, toggled with the key ``P``. The particles are
drawn at their real positions, so a particle straddling a face of the cell is
shown on one side only; this is expected, and not a sign that the periodicity is
broken.


Current status
--------------

.. warning::

   Periodic boundary conditions are implemented and usable, but they are still
   described as being in a testing phase. The work is carried out with
   *Lhassan Amarsid* and *Duc-Cuong Pham*. Check your results against a known
   case before relying on them, in particular for the deformable-cell shear.
