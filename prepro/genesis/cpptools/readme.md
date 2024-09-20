## Introduction

The program `generator` reads commands from a text-file and performs specific actions based on the commands.


## How to Use


Use the following command to run the program:

```
generator <prepro-command-file>
```
   

### Available general commands

The program currently supports the following commands:

- `open <string:filename>`
  
  Opens the file specified by `filename` (without space) for output.

- `close` 

  Closes the currently open output file, if any.

- `print <text>` 

  Prints the specified `text` to the currently open output file or to the console if no file is open.

- `addParticle <string:name> <int:group> <int:cluster> <double:homothety> <double:position_x> <double:position_y> <double:position_z> <double:angularQuaternion_s> <double:angularQuaternion_x> <double:angularQuaternion_y> <double:angularQuaternion_z>` 

  Adds a new particle with the specified parameters.


### Commands for generating a shape 

- `generateShape:sphere <string:name> <double:radius>`

  Generates a spherical shape with the specified parameters.

- `generateShape:cube <string:name> <double:radius> <vec3r:sideSize>` 

  Generates a cubic shape with the specified parameters.

- `generateShape:cuboid <name> <radius> <sideSize_x> <sideSize_y> <sideSize_z>`

  Generates a cuboid shape with the specified parameters.

- `generateShape:rectangle_xz <name> <radius> <side_x> <side_z>`

  Generates a rectangle shape (in plane xz) with the specified parameters.

- `generateShape:rhombicuboctahedron <name> <radius> <sideSize_x> <sideSize_y> <sideSize_z>`

  Generates a rhombicuboctahedron shape with the specified parameters.

- `generateShape:xyz_walls <size_x> <size_y> <size_z> <Rw>`

  Generates 3 wall-shapes with the specified parameters.


### Commands for generating a packing


- `generatePacking:wallBox <group> <LX> <LY> <LZ> <Rw>`

  Generates a packing using the wallBox method with the specified parameters.

- `generatePacking:grid <name> <origBox_x> <origBox_y> <origBox_z> <boxSize_x> <boxSize_y> <boxSize_z> <n_x> <n_y> <n_z> <group> <clust> <homothety> <randQ>`

  Generates a packing using the grid method with the specified parameters.

- `generatePacking:grid_clust <name> <origBox_x> <origBox_y> <origBox_z> <boxSize_x> <boxSize_y> <boxSize_z> <n_x> <n_y> <n_z> <group> <clust> <homothety> <randQ>`

  Generates a packing using the grid method with clustering with the specified parameters.