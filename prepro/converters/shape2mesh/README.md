# shape2mesh

Turn the r-shapes of a `.shp` file into surface meshes for 3D viewing.

Each r-shape is the Minkowski sum of its skeleton (vertices, edges, faces) with
a ball of radius `R`. `shape2mesh` meshes the exact surface of that dilated
solid — rounded edges and vertices, concavities, and open surfaces all handled —
and writes one mesh file per shape. Vertex normals are the exact analytic
normals of the surface, so even a coarse mesh shades well.

The mesh lives in the body frame of the shape (the frame stored in the `.shp`).

## Usage

```
shape2mesh <input.shp> [options]
  -f, --format obj|ply   output format            (default: obj)
  -e, --epsilon <value>  absolute sag tolerance    (default: 1% of the radius)
  -s, --simplify [deg]   merge coplanar triangles  (default angle: 0.25 degree)
  -r, --rmsh             write one <stem>.rmsh companion (all shapes in one file)
  -o, --outdir <dir>     output directory          (default: alongside input)
  -h, --help
```

## The `.rmsh` companion (viewer cache)

With `-r`/`--rmsh`, all shapes of the input are meshed into a single
`<stem>.rmsh` file (one named block per shape: positions, exact normals,
triangles). When this file sits next to a shape file, `see` draws each particle
with its skin mesh instead of the overlapping-primitive rendering. Without an
explicit `-e`, the sag is chosen per shape (0.4% of its bounding-box diagonal),
so a large flat wall and a small detailed grain both get a sensible mesh:

    shape2mesh shapes.txt --rmsh        # -> shapes.rmsh
    shape2mesh shapes.txt --rmsh -s     # smaller, flat faces merged

With `-s`, adjacent coplanar triangles are merged and the flat regions are
re-triangulated from their boundary only. The offset faces of a r-shape are
exact planes, so this removes their interior vertices without changing the
surface there: the mesh stays watertight (Euler characteristic 2, every edge on
two triangles). Curved regions (rounded edges and corners) are left untouched.
On flat-dominated shapes the triangle count drops a lot (a rounded polyhedron by
~60%, the concave L by ~50%, a plate by most of its interior); a sphere is
unchanged. The optional degree argument widens what counts as coplanar.

A region is grown by two tolerances: the normal must stay within the angle above,
and each vertex must lie within a small distance of the seed plane. That distance
is **a fraction of the requested sag** (a twentieth), not of the bounding box: it
caps how far flattening a dropped vertex can move the surface, so the merge's
residual is a fraction of the mesh's own error and shrinks with it. Measured on a
rounded polyhedron, the relative volume change stays around `1e-6` for sags from
`4e-3` down to `3e-4`. A tolerance tied to the bounding box would stay fixed while
the mesh error falls, and would climb to `2e-5` over the same range.

A file holding a single shape is written as `<inputStem>.<ext>`. A file holding
several shapes writes one file per shape, `<inputStem>_<shapeName>.<ext>`.

`epsilon` is the maximum chordal error (sag) of the mesh, in the same length
unit as the `.shp`. Smaller means finer. The grid step is derived from it as
`h ~ sqrt(8 * epsilon * R)`; halving the target sag roughly doubles the triangle
count along each dimension.

## Examples

```sh
# The demo shapes (a plate, a sphere, a rice grain)
shape2mesh ../../../examples/helloworld/shapes.shp -f obj -o /tmp

# A real polyhedron packing, as PLY, with a chosen tolerance
shape2mesh ../../../test/input/518_poly/shape.shp -f ply -e 2e-3 -o /tmp

# A Stanford bunny converted from STL
shape2mesh ../stl2shape/stanford_bunny_309_faces.shp -f obj -o /tmp
```

Open the result in Blender, MeshLab, or ParaView.

## Notes and current limits

- **Non-manifold face soups** (e.g. an imperfect STL, like the bunny) are still
  meshed as solids: the inside/outside sign falls back to ray casting, matching
  `Shape::inside()`. A genuine open surface (`isSurface` in the `.shp`) is meshed
  as a shell of thickness `2R`.
- **Large thin walls** (a big plate with a small `R`) produce very dense meshes,
  because a uniform grid fine enough for the rounded rim also tessellates the
  flat span. Coarsen with `-e` for such shapes; an analytic fast path for
  convex shapes is planned to remove this cost. A hard budget on the grid size
  keeps the tool from ever running away — beyond it the step is coarsened and a
  note is printed.
- **Detailed shapes** (a several-thousand-face STL) are handled by a bounding
  volume hierarchy for the closest-point and ray-cast queries, and the grid
  sampling runs in parallel with OpenMP. A dense STL still takes a few seconds;
  lower the resolution with `-e` for quick previews.
- The underlying mesher is validated in `test/ShapeSDF` (watertightness, exact
  on-surface vertices and normals, and second-order volume convergence).
