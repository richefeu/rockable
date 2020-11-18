## `cpptools` & the application `generator`

In the folder `prepro/genesis/cpptools`, there are some C++ pieces of code that can be used to generate shapes, packings ro entire systems. They are in the form of C++ header files, and they can be used as such.

It is possible also to use the application `generator`. This tutorial gives some examples.

## Defining a cubic shape

To define a sphero-polyhedron that is a cube, we may write the following command in a text-file named e.g. `genCube.txt`:

```
generateShape:cube Cube 0.0005 0.01
```

It means: generate a cubic shape with external size of 1 cm, and Minskowski radius of 5 mm. To generate the shape, the application `generator` is used, and from the current working directory the command is:

```
../../prepro/genesis/cpptools/generator genCube.txt > shapes.txt
```

As a result, the file `shapes.txt` should be:

```
<
name Cube
radius 0.0005
preCompDone y
nv 8
0.0045 0.0045 -0.0045
-0.0045 0.0045 -0.0045
-0.0045 -0.0045 -0.0045
0.0045 -0.0045 -0.0045
0.0045 0.0045 0.0045
-0.0045 0.0045 0.0045
-0.0045 -0.0045 0.0045
0.0045 -0.0045 0.0045
ne 12
0 1
1 2
2 3
3 0
4 5
5 6
6 7
7 4
0 4
1 5
2 6
3 7
nf 6
4 0 1 2 3
4 4 5 6 7
4 0 1 5 4
4 2 3 7 6
4 1 2 6 5
4 0 4 7 3
obb.extent 0.005 0.005 0.005
obb.e1 1 0 0
obb.e2 0 1 0
obb.e3 0 0 1
obb.center 0 0 0
volume 9.78318e-07
I/m 1.66667e-05 1.66667e-05 1.66667e-05
>
```

## Defining the wall shapes

We would like to deposite a number of cubes inside a box. So a box needs to be defined, and something similar to the cube can be done to generate a 6-wall box. Here is the file `genBox.txt`:

```
generateShape:xyz_walls   0.25 0.25 0.25   0.0005
``` 


It means: generate 3 walls that form a cubic box with inner size of 25 cm, and Minskowski radius of 5 mm. The normal to each wall is oriented towards x, y and z, and the name of the walls are respectively `x-wall`, `y-wall` and `z-wall`.  To generate the shape and to append it in the file `shapes.txt`, the command is:

```
../../prepro/genesis/cpptools/generator genBox.txt >> shapes.txt
```

## Packing the cubes inside the box

Now that the shapes are defined, the initial packing can be constructed. Here is the generation script in file `genPacking.txt`:

```
generatePacking:wallBox
1
0.25 0.25 0.25
0.0005

generatePacking:grid
Cube
0 0 0
0.25 0.25 0.25
15 15 15
0
1.0
1
```

The first command place the 6 walls that form the box, where group is 1, size is $25^3$ cm$^3$, and half-thickness is 5 mm (this latter value needs to be consistent with the one used to generate the wall shapes).

The second command place $15^3$ particles of shape named `Cube` inside a box ($25\times25\times25$ cm$^2$), from origin point (0, 0, 0). The group of the added particles is 0, the size factor (homothety) is 1.0, and the last value says that the cubes will have random orientations.

The generator is invoked that way:

```
../../prepro/genesis/cpptools/generator genPacking.txt > sample.txt
```

And we get in the terminal the number of added bodies for each command:

```
@wallBox, Number of added bodies: 6
@grid, Number of added bodies: 3375
```

## Assembling the input conf-file

Now, we have all the pieces to assemble the input file. Let's start with the file `input_head.txt`: `cp input_head.txt input.txt`.

Open this file and edit the following lines:

```
shapeFile shapes.txt
nDriven 6
Particles 3381
```

Notice that the number of particles is the total number 3375 cubes + 6 walls (placed at the beginning of the list).
The content of the file `sample.txt` can be copied at the end of the file `input.txt`.

The simulation is ready to be run:

```
../../src/run input.txt
```