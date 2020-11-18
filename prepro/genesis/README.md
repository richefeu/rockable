This folder holds tools to generate shapes and/or packings.

## File names

- Each file named `cpptools/generateShape_*.hpp` hold a function (potentially several functions) that will generate one or more shape(s) in the shape-format of `Rockable`.

- Each file named `cpptools/generatePacking_*.hpp` is useful to generate packing as the name suggests. It often uses some shape-generator to produce a list of particles to be included in the input files of `Rockable`.

- the file `cpptools/addParticle.hpp` provides a function useful to add a particle (in the right format) by giving it a shape-reference, a group identifier, a cluster identifier, an homothety, a position and an orientation (quaternion).

- The folders (*e.g.*, `SpherePacker`) are some system generators. They are applications that use the shape and packing generators to build system ready for being used by `Rockable`.

## Shape generators

- `sphere`: a sphero-polyhedron with a single node 
- `cube`: a cube with a given external side size
- `cylindricalMold`: cylindrical container with axis aligned with the $y$ direction 
- `xyz_walls`: 3 axis-aligned walls with chosen thickness
- `stick`: a sphero-polyhedron with 2 nodes connected with a single edge along the $x$ direction


## Packing generators

TODO

## System generators (shapes & packing)

- `SpherePacker`: it packs, with a targetted solid fraction, spherical particles inside a container (6-walls or cylinder)
- `TubePacker`: it packs, with a targetted solid fraction, sphero-line particles inside a container (6-walls or cylinder)