The source code can be cloned from its repository (provided that you've got the right)

```
git clone https://<sourceSupLogin>@git.renater.fr/authscm/<sourceSupLogin>/git/rockable/rockable.git
```

If you want to compile with the `makefile`, you will also need to clone toofus into your home directory:

```
cd ~
git clone https://www.github.com/richefeu/toofus.git
```

For compilations, `gcc` (version >=10 if possible) is required.
Some dependencies need to be installed (the examples are for a linux system):

- for `run`(or `rockable`):
 
  - `spdlog`
  
  ```
  sudo apt-get install libspdlog-dev
  ```
  
  - `tclap`
  
  ```
  sudo apt-get install libtclap-dev
  ```

- for `see`:


  - `freeglut3-dev`

  ```
  sudo apt-get install freeglut3-dev 
  ```

  - `libxmu-dev`
  
  ```
  sudo apt-get install libxmu-dev
  ```

  - `libxi-dev`

  ```
  sudo apt-get install libxi-dev
  ```
  
- for `see3` (`visuedit`):

  - `libglfw3-dev`  

  ```
  sudo apt-get install libglfw3-dev
  ```

### How to compile (with `make`)

```
make -j 8
```

### How to compile (with `cmake`)
 
```
mkdir build
cd build
cmake ..
make -j 8
```

#### Disable `see`

```
cmake .. -DROCKABLE_COMPILE_SEE=OFF
```

#### How to install with profiling tools

```
cmake .. -DROCKABLE_ENABLE_PROFILING=ON
```
