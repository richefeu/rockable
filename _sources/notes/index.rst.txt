.. _technicalNotes:

Technical notes
===============

Some parts of Rockable rest on methods that need more than a reference page:
derivations, choices made in the implementation, validation and known limits.
They are documented in technical notes, written like short articles, with
their references to the literature. Each note is a PDF built from the LaTeX
sources kept in this folder of the repository
(``sphinxdoc/source/notes/<note>/``, built with ``make``).

.. contents::
   :local:


.. _noteSkinMeshes:

Watertight skin meshes of three-dimensional R-shapes
----------------------------------------------------

*Exact distance, visualisation and mass properties* --- V. Richefeu, October 2026.

:download:`Download the note (PDF) <skin-meshes.pdf>`

A sphero-polyhedron is stored as a skeleton and a radius, so its surface is
never built, and the viewers draw it with overlapping spheres, tubes and plates,
which fails with transparency and cannot be exported. The note documents the
tools that build the surface explicitly:

- the exact signed distance to an R-shape (the distance to the skeleton minus
  the radius) and its analytic gradient, with the three rules that give its sign
  (clean manifold, imported face soup, open surface) --- class ``ShapeSDF``;
- the watertight skin mesh, by marching cubes and Newton projection onto the
  exact surface, driven by a single sag tolerance --- ``buildShapeMesh``;
- the coplanar merge that removes the extra triangles of the flat faces without
  moving the surface --- ``simplifyCoplanar``;
- the volume, centre of mass and inertia by the divergence theorem on the mesh,
  compared with the Monte-Carlo estimate of ``Shape::massProperties``;
- the command-line mesher ``shape2mesh`` and the ``.rmsh`` file read by ``see``
  and ``shapeSurvey`` (see also :ref:`converters`), the validation tests in
  ``test/ShapeSDF``, and stress tests on difficult skeletons, with the cases
  where the mesher fails.

.. note::

   In the library, the mesh gives the volume and the area of a shape. The centre
   of mass and the inertia computed from the mesh are described in the note but
   are not yet used by ``Shape::massProperties``, which still relies on the
   Monte-Carlo integration (keyword ``MCnstep`` of the shape file).
