`stl2Shape` allows for the conversion from a **binary** STL file to a SHAPE file that can be used by `Rockable`.
It uses the library `tclap` to parse the command line arguments. You can thus use the `-h` or `--help` argument to show all possible arguments: 

```
USAGE: 

   ./stl2Shape  [-m <double>] [-s <double>] -r <double> -i <string> [--]
                [--version] [-h]


Where: 

   -m <double>,  --maxLength <double>
     Set the maximum length of the output object

   -s <double>,  --scaleFactor <double>
     Set the re-scale factor of the output object

   -r <double>,  --radius <double>
     (required)  Minskowski radius to be used

   -i <string>,  --input <string>
     (required)  Name to STL binary file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Convert a binary STL file to a shape that can be used by Rockable
```

The arguments `--maxLength` and `--scaleFactor` will both affect the chosen radius. No extra-data (OBB, inertia, volume, etc.) will be saved in the output shape-file. this means that the application `shapeSurvey` will be used to compute these extra-data.
