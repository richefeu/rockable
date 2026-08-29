#import "common.typ": *

#show: rockable-sheet.with(
  title: [`sweepable`, parameter sweeps and batch runs],
)

= Principle

#concept-block(body: [
  `sweepable` reads a script and generates a *tree of input files*, one per
  parameter combination, by editing a template line by line. It knows nothing
  about `Rockable`: it manipulates text files, so it works with any code whose
  input is a text file.

  Two pieces, no external dependency, C++17:
  #kwtable(
    file("simuFileManips.hpp"), [The header-only engine (class `SFManip`).],
    file("sweepable.cpp"), [The script interpreter. Build it with `make`.],
    file("batchable.sh"), [The companion runner that launches the generated
      inputs in parallel.],
  )
])

= The three steps

#concept-block(body: [
  #inline("1. Build")
  ```sh
  make
  ```

  #inline("2. Generate the runner")
  Running `sweepable` *without argument* writes `batchable.sh` into the current
  folder:
  ```sh
  ./sweepable
  chmod +x batchable.sh
  ```

  #inline("3. Generate the inputs, then run them")
  ```sh
  ./sweepable sweep.txt
  ./batchable.sh generated.txt "rockable -j 2" 8
  ```
])

= `batchable.sh`

#concept-block(body: [
  #align(center, command("batchable.sh <input_list> <command> <nprocs>"))

  #kwtable(
    command("input_list"), [The file written by `.saveCollection`: one generated
      input path per line.],
    command("command"), [The executable to run. It receives the *file name* as
      its argument.],
    command("nprocs"), [Maximum number of jobs run in parallel.],
  )

  Each job `cd`s into the folder of its input file, runs the command there, and
  redirects everything into a #file("run.log") of that folder. So the conf-files
  of each case stay in their own directory.

  #note[Count the cores: `rockable -j 2` times 8 parallel jobs already asks for
  16 threads.]
])

= Script syntax

#concept-block(body: [
  Commands start with a dot, arguments are space-separated, quotes preserve
  spaces, and `#` starts a comment.

  ```txt
  # This is a comment
  .read input.txt
  .appendLine "# generated automatically"
  ```

  #inline("Variables")
  ```txt
  .set dt 0.5
  .createFolder "output/dt_$dt"
  .createFolder "output/${dt}_run"
  ```
  Use `$name`, or `${name}` when the name is followed by other characters.

  #inline("Loops")
  ```txt
  .foreach variable value1 value2 ...
      ...
  .endforeach
  ```
  Loops nest freely. The loop variable is restored to its previous value when
  the loop exits.

  #inline("Verbosity")
  #command(".verbose") (default), #command(".mute") or #command(".quiet").
])

= Command reference

#concept-block(body: [
  #inline("File I/O")
  - #command(".read <file>") ~Load a file into memory.
  - #command(".setOutputFilename <name>") ~Change the output name (default
    `input.txt`).
  - #command(".saveInFolder") ~Save into the folder remembered by
    `.createFolder`.
  - #command(".saveInFolder \"<folder>\"") ~Save into an explicit folder,
    creating it if needed.

  #inline("Selecting and editing lines")
  - #command(".findLineStartingWith <prefix>") ~Select the first matching line.
  - #command(".replaceBy \"<text>\"") ~Replace the selected line.
  - #command(".replaceAllStartingWith <prefix> \"<text>\"") ~Replace every
    matching line.
  - #command(".replaceInLine <old> <new>") ~Replace a substring inside the
    selected line.
  - #command(".insertAfter <prefix> \"<line>\"") ~Insert a line after the first
    match.
  - #command(".appendLine \"<line>\"") ~Append a line at the end.

  #inline("Folders and collection")
  - #command(".createFolder \"<path>\"") ~Create a directory and remember it.
  - #command(".clearCollection") ~Empty the list of generated files.
  - #command(".saveCollection <file>") ~Write that list to disk --- this is what
    `batchable.sh` consumes.

  #inline("Error handling")
  - #command(".silentIfNotFound") ~Ignore a missing match instead of throwing.
  - #command(".silentIfNotFound false") ~Restore the errors.
])

= A complete sweep

#concept-block(body: [
  #file("sweep.txt"):
  ```txt
  .foreach nd 0.0 0.02 0.05
      .foreach dt 1e-6 5e-7

          .read input.txt

          .findLineStartingWith numericalDampingCoeff
          .replaceBy "numericalDampingCoeff $nd"
          .findLineStartingWith dt
          .replaceBy "dt $dt"

          .createFolder "runs/nd_${nd}_dt_$dt"
          .saveInFolder

      .endforeach
  .endforeach

  .saveCollection generated.txt
  ```

  produces
  ```txt
  runs/nd_0.0_dt_1e-6/input.txt
  runs/nd_0.0_dt_5e-7/input.txt
  runs/nd_0.02_dt_1e-6/input.txt
  ...
  ```

  #note[`.findLineStartingWith` matches the *first* line only. For a keyword
  that appears once per group pair --- `muContact 0 1`, `muContact 1 1` --- match
  the full prefix with `.replaceAllStartingWith "muContact 0 1"`, or the other
  pairs will be overwritten too.]

  #note[The shape file and the companion files are *not* copied. Either copy
  them into each folder yourself, or leave the template pointing at a shared
  path.]
])

= The C++ API

#concept-block(body: [
  `SFManip` can be used directly, without the interpreter. Every method returns
  `SFManip&`, so the calls chain, and the "replace" family is printf-style.

  ```cpp
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
  ```

  #inline("Methods")
  `read(filename)`, `setOutputFilename(name)`,
  `findLineStartingWith(prefix)`, `replaceBy(fmt, ...)`,
  `replaceAllStartingWith(prefix, fmt, ...)`,
  `replaceInLine(target, replacement)`, `insertAfter(prefix, newLine)`,
  `appendLine(line)`, `createFolder(fmt, ...)`, `saveInFolder()`,
  `saveInFolder(fmt, ...)`, `silentIfNotFound(bool)`.

  #inline("Static collection")
  ```cpp
  SFManip::Collection         // std::vector<std::string> of paths
  SFManip::clearCollection()
  SFManip::saveCollection(file)
  ```
])

= A `Rockable` study, end to end

#concept-block(body: [
  ```sh
  # 1. one template folder
  mkdir template && cd template
  #    input.txt, shapes.txt, drivingSystem.txt ...

  # 2. generate the cases
  cd ..
  ./sweepable sweep.txt        # -> runs/*/input.txt + generated.txt

  # 3. give each case what it needs beside its input
  while read f; do
    cp template/shapes.txt template/drivingSystem.txt "$(dirname "$f")"
  done < generated.txt

  # 4. run them, 8 at a time
  ./batchable.sh generated.txt "rockable -j 2 -v 3" 8

  # 5. collect one extractor output per case
  for d in runs/*/; do
    echo "$d $(tail -1 "$d/kineticEnergy.txt")"
  done
  ```

  #inline("Why a template folder")
  `sweepable` only writes the file it read. Everything else a simulation needs
  --- the shape library, #file("drivingSystem.txt"),
  #file("dataExtractors.txt") --- has to be placed next to each generated input,
  which is what step 3 does.

  #inline("Checking before launching")
  ```sh
  wc -l generated.txt          # how many cases
  head -3 generated.txt
  diff template/input.txt runs/nd_0.0_dt_1e-6/input.txt
  ```
  The `diff` is the quickest way to confirm that the intended lines --- and only
  those --- were changed.
])
