.. _softParticles:

Deformable particles
====================

By default a Rockable particle is rigid: it has six degrees of freedom, three of
translation and three of rotation. The *soft particles* extension adds a
**homogeneous deformation** to each body, so that a particle can change shape
under load without being discretised into a mesh.

The intent, stated in the design of the feature, is not to resolve a full
kinematic field with many degrees of freedom, but to make a few of them
possible: enough to model soft grains, polymers or biological tissues, while
keeping the cost of a discrete element method.

.. important::

   This feature requires

   .. code-block:: sh

      cmake .. -DROCKABLE_ENABLE_SOFT_PARTICLES=ON

   Without it, the keyword is accepted but only prints
   ``SOFT_PARTICLES_NOT_ENABLED when Rockable was compiled``.


Activating the feature
----------------------

- ``useSoftParticles`` (*double*) **Young** (*double*) **Poisson**

  Activates the homogeneous straining of the particles and builds the
  compliance from the Young modulus :math:`E` and the Poisson ratio
  :math:`\nu` of the material the bodies are made of.

.. code-block:: text
   :caption: input.txt

   useSoftParticles 1e6 0.0

Both parameters are global: all the bodies share the same compliance, whatever
their group. A Poisson ratio of zero, as above, means the body shortens along
the loading direction without bulging sideways, which is the simplest case to
interpret.

.. note::

   The compliance relates the average stress carried by a body to its uniform
   strain. The stiffness of the **contacts** is a different matter, set by
   ``knContact`` and ``ktContact`` as usual: a soft particle in a stiff contact
   and a stiff particle in a soft contact are two different models, and both
   are available.


Current status
--------------

.. warning::

   The extension is under active development, in the postdoctoral work of
   *Mukesh Singh Bisht*, and it is not yet complete. In particular, the uniform
   transformation of the particles is **not written to the conf-files**: the
   code that saved it is commented out, and the corresponding field is not read
   back either.

   The practical consequence is that a computation using soft particles cannot
   be restarted exactly from a dump. The positions, velocities and orientations
   are restored, but the accumulated deformation of each body is lost, and the
   bodies restart undeformed. Run such a simulation in one go, or check the
   state of the feature in the source before relying on a restart.

.. tip::

   ``rockable -b`` prints ``ROCKABLE_ENABLE_SOFT_PARTICLES`` among the
   compilation options, which is the quickest way to tell whether the binary at
   hand supports the feature.
