`stl2Shape` allows for the conversion from a **binary** STL file to a SHAPE file that can be used by `Rockable`.
It uses the library `tclap` to parse the command line arguments. You can thus use the `-h` or `--help` argument to show all possible arguments: 

```

USAGE: 

   stl2shape  [-chz] [--version] [-m <double>] [-n <string>] [-o <string>]
              [-s <double>] -r <double> input


Where: 

   -r <double>,  --radius <double>
     (required) Radius of the rounded edges (Minkowski radius) to be used

   -s <double>,  --scaleFactor <double>
     Set the (re-)scale factor of the output object

   -m <double>,  --maxLength <double>
     Set the maximum length (sieving size) of the output object

   -c,  --clean
     Remove duplicated edges

   -z,  --scaleRadius
     Rescale the radius

   -n <string>,  --name <string>
     Name of the shape (default: the name of the STL file, without its
     folder)

   -o <string>,  --output <string>
     Name of the output shape file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   input <stl file>
     (required) Name of the input STL binary file

   Convert a binary STL file to a shape that can be used by Rockable

```

The options `-n` and `-o` give the name of the shape (by default, the name of the STL file, without its folder) and the name of the file written (by default `mesh.shp`). The arguments `-m, --maxLength` and `-s, --scaleFactor` will both affect the chosen radius if option `-z, --scaleRadius` is used. No extra-data (OBB, inertia, volume, etc.) will be saved in the output shape-file. this means that the application `shapeSurvey` will be used to compute these extra-data.
