## Introduction

The program `generator` reads commands from a text-file and performs specific actions based on the commands.


## How to Use


Use the following command to run the program:

```
generator <prepro-command-file>
```
   

### Available general commands

The program currently supports the following commands:

- `open [string:filename]`
  
  Opens the file specified by `filename` (without space) for output.

- `close` 

  Closes the currently open output file, if any.

- `print [text up to the line end]` 

  Prints the specified `text` to the currently open output file or to the console if no file is open. The spaces at the beginning and at the end of the text are trimmed. After printing, a cariage return is added (next print will star at the next line).
  
- `print> [text up to the line end]` 

  Similar to `print` but without cariage return at the end. Next Print will follow at the end.

- `<compute`

  Blabla
  
- `compute`

- `transformation:rotate`

  Blabla

- `transformation:translate`
  
  Blabla

- `transformation:reset`
  
  Blabla
  
- `addParticle [string:name] [int:group] [int:cluster] [double:homothety] [double:position_x] [double:position_y] [double:position_z] [double:angularQuaternion_s] [double:angularQuaternion_x] [double:angularQuaternion_y] [double:angularQuaternion_z]` 

  Adds a new particle with the specified parameters.


### Commands for generating a shape 

- `generateShape:sphere [string:name] [double:radius]`

  Generates a spherical shape with the specified parameters.

- `generateShape:cube [string:name] [double:radius] [double:sideSize_x] [double:sideSize_y] [double:sideSize_z]` 

  Generates a cubic shape with the specified parameters.

- `generateShape:cuboid [string:name] [double:radius] [double:sideSize_x] [double:sideSize_y] [double:sideSize_z]`

  Generates a cuboid shape with the specified parameters.

- `generateShape:rectangle_xz [string:name] [double:radius] [double:side_x] [double:side_z]`

  Generates a rectangle shape (in plane xz) with the specified parameters.

- `generateShape:rhombicuboctahedron [string:name] [double:radius] [double:sideSize_x] [double:sideSize_y] [double:sideSize_z]`

  Generates a rhombicuboctahedron shape with the specified parameters.

- `generateShape:xyz_walls [double:size_x] [double:size_y] [double:size_z] [double:Rw]`

  Generates 3 wall-shapes with the specified parameters.


### Commands for generating a packing


- `generatePacking:wallBox [group] [double:LX] [double:LY] [double:LZ] [double:Rw]`

  Generates a packing using the wallBox method with the specified parameters.

- `generatePacking:grid [string:name] [double:origBox_x] [double:origBox_y] [double:origBox_z] [double:boxSize_x] [double:boxSize_y] [double:boxSize_z] [int:n_x] [int:n_y] [int:n_z] [int:group] [int:clust] [double:homothety] [int:randQ]`

  Generates a packing using the grid method with the specified parameters.

- `generatePacking:grid_clust [string:name] [double:origBox_x] [double:origBox_y] [double:origBox_z] [double:boxSize_x] [double:boxSize_y] [double:boxSize_z] [int:n_x] [int:n_y] [int:n_z] [int:group] [int:clust] [double:homothety] [int:randQ]`

  Generates a packing using the grid method with clustering with the specified parameters.