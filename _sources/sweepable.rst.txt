Application `sweepable`
=======================

``sweepable`` is a lightweight command-line tool for simulation parameter sweeps.
It reads a script file and generates a tree of input files — one per parameter
combination — ready to be launched in parallel with the companion script
``batchable.sh``.

Both tools ship as a single self-contained file:

- ``simuFileManips.hpp`` — header-only C++ library (the engine)
- ``sweepable.cpp`` — the script interpreter (build with ``make``)

No external dependencies. Requires C++17.

.. contents::
   :local:
   :depth: 2


Quick Start
-----------

Build
^^^^^

.. code-block:: sh

   make

Generate the batch runner
^^^^^^^^^^^^^^^^^^^^^^^^^

Running ``sweepable`` without arguments writes ``batchable.sh`` to the current
directory:

.. code-block:: sh

   ./sweepable
   # → File 'batchable.sh' has been created
   chmod +x batchable.sh

Write a sweep script
^^^^^^^^^^^^^^^^^^^^

``sweep.txt``:

.. code-block:: text

   .foreach solver cg gmres
       .foreach dt 0.1 0.5 1.0

           .read input.txt

           .findLineStartingWith solver
           .replaceBy "solver $solver"

           .findLineStartingWith dt
           .replaceBy "dt $dt"

           .createFolder "output/${solver}_dt_$dt"
           .saveInFolder

       .endforeach
   .endforeach

   .saveCollection generated.txt

Run the sweep
^^^^^^^^^^^^^

.. code-block:: sh

   ./sweepable sweep.txt

This produces:

.. code-block:: text

   output/cg_dt_0.1/input.txt
   output/cg_dt_0.5/input.txt
   output/cg_dt_1.0/input.txt
   output/gmres_dt_0.1/input.txt
   output/gmres_dt_0.5/input.txt
   output/gmres_dt_1.0/input.txt
   generated.txt

Launch all runs in parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: sh

   ./batchable.sh generated.txt mysim 4

``batchable.sh`` reads ``generated.txt`` (produced by ``.saveCollection``),
runs ``mysim input.txt`` inside each folder, and uses up to 4 parallel
processes. Output from each run is captured in ``run.log`` inside its folder.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Argument
     - Description
   * - ``input_list``
     - File produced by ``.saveCollection``
   * - ``command``
     - Executable to run (receives the input filename as argument)
   * - ``nprocs``
     - Maximum number of parallel jobs


Script Language
---------------

General syntax
^^^^^^^^^^^^^^

Commands start with a dot. Arguments are space-separated. Use quotes to
preserve spaces. Comments start with ``#``.

.. code-block:: text

   # This is a comment
   .read input.txt
   .appendLine "# generated automatically"

Variables
^^^^^^^^^

Define a variable with ``.set``:

.. code-block:: text

   .set name value

Use it with ``$name`` or ``${name}`` (braces avoid ambiguity in composite strings):

.. code-block:: text

   .set dt 0.5
   .createFolder "output/dt_$dt"
   .createFolder "output/${dt}_run"

Loops
^^^^^

.. code-block:: text

   .foreach variable value1 value2 ...
       ...
   .endforeach

Loops are fully nestable. The loop variable is restored to its previous value
after the loop exits.

Verbosity
^^^^^^^^^

.. code-block:: text

   .verbose    # enable output (default)
   .mute       # suppress output
   .quiet      # same as .mute


Command Reference
-----------------

File I/O
^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Command
     - Description
   * - ``.read <file>``
     - Load a file into memory
   * - ``.setOutputFilename <name>``
     - Change output filename (default: ``input.txt``)
   * - ``.saveInFolder``
     - Save to the current folder
   * - ``.saveInFolder "<folder>"``
     - Save to a specific folder (creates it if needed)

Line Selection and Editing
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 52 48

   * - Command
     - Description
   * - ``.findLineStartingWith <prefix>``
     - Select the first matching line
   * - ``.replaceBy "<text>"``
     - Replace the selected line
   * - ``.replaceAllStartingWith <prefix> "<text>"``
     - Replace all matching lines
   * - ``.replaceInLine <old> <new>``
     - Replace a substring in the selected line
   * - ``.insertAfter <prefix> "<line>"``
     - Insert a line after the first match
   * - ``.appendLine "<line>"``
     - Append a line at the end of the file

Folders
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Command
     - Description
   * - ``.createFolder "<path>"``
     - Create a directory and remember it for ``.saveInFolder``

Collection
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Command
     - Description
   * - ``.clearCollection``
     - Clear the list of generated files
   * - ``.saveCollection <file>``
     - Write the list of generated files to disk

Error Handling
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Command
     - Description
   * - ``.silentIfNotFound``
     - Ignore missing matches instead of throwing
   * - ``.silentIfNotFound false``
     - Re-enable errors for missing matches


C++ API
-------

``SFManip`` can also be used directly in C++ without the script interpreter.

Basic example
^^^^^^^^^^^^^

.. code-block:: cpp

   #include "simuFileManips.hpp"

   int main() {
       SFManip sf;
       sf.read("input.txt")
         .findLineStartingWith("dt")
         .replaceBy("dt %.3f", 0.5)
         .createFolder("output/run1")
         .saveInFolder();
   }

Full sweep example
^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   #include "simuFileManips.hpp"

   int main() {
       SFManip::clearCollection();

       for (const auto& solver : {"cg", "gmres"}) {
           for (double dt : {0.1, 0.5, 1.0}) {
               SFManip()
                 .read("input.txt")
                 .findLineStartingWith("solver")
                 .replaceBy("solver %s", solver)
                 .findLineStartingWith("dt")
                 .replaceBy("dt %.3f", dt)
                 .createFolder("output/%s_dt_%.3f", solver, dt)
                 .saveInFolder();
           }
       }

       SFManip::saveCollection("generated.txt");
   }

Method summary
^^^^^^^^^^^^^^

All methods return ``SFManip&`` for chaining.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``read(filename)``
     - Load file into memory
   * - ``setOutputFilename(name)``
     - Set output filename
   * - ``findLineStartingWith(prefix)``
     - Select first matching line
   * - ``replaceBy(fmt, ...)``
     - Replace selected line (printf-style)
   * - ``replaceAllStartingWith(prefix, fmt, ...)``
     - Replace all matching lines
   * - ``replaceInLine(target, replacement)``
     - Replace substring in selected line
   * - ``insertAfter(prefix, newLine)``
     - Insert line after match
   * - ``appendLine(line)``
     - Append line at end
   * - ``createFolder(fmt, ...)``
     - Create directory and store it
   * - ``saveInFolder()``
     - Save to stored directory
   * - ``saveInFolder(fmt, ...)``
     - Save to explicit directory
   * - ``silentIfNotFound(bool)``
     - Suppress errors for missing matches

Static collection API
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   SFManip::Collection            // std::vector<std::string> of generated paths
   SFManip::clearCollection()     // clear the list
   SFManip::saveCollection(path)  // write list to file