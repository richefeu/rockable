## Introduction

The program `generator` reads commands from a text-file and performs specific actions based on the commands.


## How to Use


Use the following command to run the program:

```
generator <prepro-command-file>
```

As an example, here is a command script:

```
# -------------------
open shapes.txt
generateShape:cuboid Cuboid 0.0005 0.03 0.01 0.01
generateShape:xyz_walls   0.25 0.25 0.25   0.0005
close

# -------------------
open input.txt
print Rockable 21-08-2022
print #### time #############################
print t 0
print tmax 2.0
print dt 5e-6
print interVerlet 0.001
print interConf 0.01
print 
print #### options ##########################
print forceLaw Default
print Integrator Beeman
print 
print #### proximity of the particles #######
print DVerlet 0.005
print dVerlet 0.0025
print 
print #### properties for particles #########
print density 0 2700
print density 1 2700
print 
print #### interaction parameters ###########
print knContact 0 1 1e6
print en2Contact 0 1 0.02
print ktContact 0 1 1e6
print muContact 0 1 0.6
print krContact 0 1 0.0
print murContact 0 1 0.0
print 
print knContact 0 0 1e6
print en2Contact 0 0 0.02
print ktContact 0 0 1e6
print muContact 0 0 0.6
print krContact 0 0 0.0
print murContact 0 0 0.0
print 
print #### the system #######################
print iconf 0
print shapeFile shapes.txt
print nDriven 0
#incrementNoParticles 1
Particles:placeHolder


generatePacking:RandomClosePacking
Cuboid
CUBOID Y
0.2 0.2 0.3
0.03 0.01 0.01
1.0 1.0
2000
0.9
0
0

close
```


Remember to leave at least one blank line at the end of the file to avoid issues with interpreting the command file. This issue may be fixed someday :)


### Available general commands

The program currently supports the following commands:

- `open [string:filename]`
  
  Opens the file specified by `filename` (without space) for output. If no file is open, the output will be directed to the console.

- `close` 

  Closes the currently open output file, if any.

- `print [text up to the line end]` 

  Prints the specified `text` to the currently open output file or to the console if no file is open. The spaces at the beginning and at the end of the text are trimmed. After printing, a cariage return is added (next print will star at the next line).
  
  This feature is particularly useful for printing raw data into shape-files or configuration files.
  
- `print> [text up to the line end]` 

  Similar to `print`, but without adding a carriage return at the end. The next print statement will continue on the same line.


  
- `compute [text before the result] [operation]` 

  This command functions similarly to the print command but includes the ability to parse a simple mathematical expression. Note that the parser is quite basic.

  Example of use: `compute Particles 5^3 + 6` will ouput `Particles 131`
   
- `<compute [operation]`

  Similar to `compute` but without text before the operation. Example of use:
  
  ```
  print> knContact 0 0
  <compute 500000 * 2
  ```
  
  Thiswill generate `knContact 0 0 1000000`

- `Particles:placeHolder`

  When this command is added, the keyword `Particles` followed by the tracked number of particles will be printed when the file is closed. The particle count is updated each time a command that adds particles is used. Additionally, the number of particles can be increased using the command `incrementNoParticles.

- `incrementNoParticles [number]`

  Increases the tracked number of particles by `number` (integer). This command is useful for manually updating the particle count, which will be printed when the file is closed if the `Particles:placeHolder command is used.



When executing the script sequentially, a global transformation is maintained and utilized by most of the packing commands. The transformation manipulations are performed by the following commands:

- `transformation:add_translation`
  
  This command adds a translation component to the current transformation, allowing you to shift the position of the object in space.

- `transformation:add_rotation`

  This command adds a rotation component to the current transformation, enabling you to rotate the object around a specified axis.

- `transformation:reset`
  
  This command resets the transformation to its initial state, effectively removing any translations or rotations that have been applied. It is also possible to reset only translation or rotation using respectively `tranformation:reset_translation` or `tranformation:reset_rotation`


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


- `addParticle [string:name] [int:group] [int:cluster] [double:homothety] [double:position_x] [double:position_y] [double:position_z] [double:angularQuaternion_s] [double:angularQuaternion_x] [double:angularQuaternion_y] [double:angularQuaternion_z]` 

  Adds a new particle with the specified parameters.


- `generatePacking:wallBox [group] [double:LX] [double:LY] [double:LZ] [double:Rw]`

  Generates a packing using the wallBox method with the specified parameters.

- `generatePacking:grid [string:name] [double:origBox_x] [double:origBox_y] [double:origBox_z] [double:boxSize_x] [double:boxSize_y] [double:boxSize_z] [int:n_x] [int:n_y] [int:n_z] [int:group] [int:clust] [double:homothety] [int:randQ]`

  Generates a packing using the grid method with the specified parameters.

- `generatePacking:grid_clust [string:name] [double:origBox_x] [double:origBox_y] [double:origBox_z] [double:boxSize_x] [double:boxSize_y] [double:boxSize_z] [int:n_x] [int:n_y] [int:n_z] [int:group] [int:clust] [double:homothety] [int:randQ]`

  Generates a packing using the grid method with clustering with the specified parameters.