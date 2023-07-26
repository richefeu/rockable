## What is `Rockable`?

Rockable is a DEM code written in C++ by <vincent.richefeu@3sr-grenoble.fr>. The two main specificities of the code are (_i_) to hold sphero-polyhedral shapes, (_ii_) to manage breakable interfaces. It is developed for an **academic usage**. This means that the code is not intended to be a tool for all purposes. It can easily be used to do what it is designed for, but to extend it, it is necessary to master both the model (DEM, complex shapes and interaction laws) and its implementation (data structure). The benefit of a good understanding is to avoid a "hacking" that would eventually limit the developed possibilities. In other words, the design of the code (neither too specific nor too general) is intended to avoid any tendency towards a single thought.

The program uses very few (if any) third-party libraries. The calculation code itself is self-sufficient. It relies only on a few small "headers-only" libraries embedded in the folder `common`. On the other hand, the 3D visualization software (`see`) uses the OpenGL programming interface.

The use of the code is not interfaced by any tool (like lua, python or any graphical interface) to facilitate its use, except the input format as described [here](https://richefeu.gitbook.io/cdm/dem/format-of-configuration-files-conf-files). This makes it particularly streamlined and greatly facilitates its integration with other calculation codes. It is in this sense that Rockable is qualified of "academic code".

## Folders

* `common`: header-only libraries used by Rockable or processing tools. This is now replace by the library `toofus` that can be found on [GitHub](https://richefeu.github.io/toofus/).
* `doc`: place to generate the doxygen html documentation of the source files of Rockable
* `sphinxdoc`: user documentation (sphinx with ReStructuredText)
* `examples`: examples for usage tutorials or for testing features
* `prepro`: some pre-processing tools
* `src`: Rockable C++ source files

## Credits

The code was initially developed by _Vincent Richefeu_, at Laboratoire 3SR, to model rockfalls and rock avalanches. This has been done through The PhD work of _Stiven Cuervo_ and _Bruna Garcia_, but it actually started before in a code named DEMbox (no longer maintained).

Then, the breakable interfaces have been implemented during the PhD work of _Marta Stasiak_. A number of improvements have been added at that time thanks to intensive review with _Gael Combe_, Laboratoire 3SR.

New functionalities are being studied thanks to new collaborations of people from CEA, IATE and CNRS. For example, _Lhassan Amarsid_ (CEA) is working on the introduction of periodic boundary conditions, and multi-processor computing with domain decomposition. _Farhang Radjai_ and students, may introduce new breakable interfaces with energy-based criteria. 

## How to install

The source code can be cloned from the SourceSup git repository. Obvioulsly, you need to be allowed. The git command may look something like:

```sh
git clone https://yourLogin@git.renater.fr/authscm/yourLogin/git/rockable/rockable.git
```

where `yourLogin` is your own login in SourceSup. Another way to get this command line is to go in "Ma Page" in SourceSup site, find "Rockable" and then "Code Source". The command can be copy-pasted (choose "via smart HTTP").

`Rockable` uses a header-only library named `toofus`. It needs to be downloaded from [https://richefeu.github.io/toofus/](https://richefeu.github.io/toofus/) (ZIP, TAR or *via* GitHub), and unzipped somewhere. The installation consists simply to copy all the header files in `/usr/local/include/toofus`. To do so, type `make` in `toofus`' folder.

Other libraries are required and need to be installed using the OS package manager (use `brew` for apple). They are `tclap`, `glfw3`, `libpng2` (optionnal).

If you are lucky, the compilation is as simple as `make` in the `src`folder of `Rockable` sources.


## Some envisioned features

- should have a look at [tinyexpr](https://github.com/codeplea/tinyexpr) so that the user can provide some custom expressions
- processing application could be written using [libcmdf](https://github.com/ronen25/libcmdf)
- https://github.com/p-ranav/indicators