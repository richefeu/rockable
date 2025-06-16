### Installation of `neper`

The sources need first to be downloaded from the site: [https://neper.info]()

To quickly compile and install it, the following solutions may be used (example for MacOS):

1. the library `gsl` needs to be installed: `brew install gsl`
2. install `nlopt` is required for graingrowth morpho: `brew install nlopt`
3. since the default compiler is `clang`, and openMP is not okay, we need to install gcc (eg. version 10)

Here are the step-by-step commands

```sh
cd src
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++-10 -DCMAKE_C_COMPILER=gcc-10
make
make install
``` 




### Using `neper` to build the voronoi polyhedra



```sh
neper -T -n 100 -domain "sphere(1,1000)" -morpho graingrowth 
neper -V n100-id1.tess -datacellcol id -datacelltrs 0.5 -cameraangle 12 -imagesize 600:600 -showedge "polynb>1" -cameraangle 9 -print preview
```


```
neper -T -n 100 -domain "sphere(1,1000)" -format stl:bycell
```