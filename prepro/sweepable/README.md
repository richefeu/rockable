
# C++17 class SFManip and command processor `sweepable`

SFManip is a lightweight C++ text-file manipulation library together with a command processor application (it uses a small command-line language) designed for simulation preprocessing, parameter sweeps, and automated generation of input files.

The project provides:

- A fluent C++ API (`SFManip`)
- A script interpreter `sweepable` for batch workflows with:
  - Variable substitution
  - Nested loops
  - Automatic folder creation
  - Generated-file collection tracking

The system is especially useful for sensibility studies.



## Design Philosophy

SFManip is intentionally:

* lightweight
* dependency-free
* line-oriented
* simulation-oriented
* easy to use

It is not intended to replace:

* full scripting languages
* regex engines
* parsers
* template systems

Instead, it provides a compact workflow tool for automated generation of structured text inputs.



## Features

### C++ API

- Read and modify text files
- Replace lines
- Insert lines
- Replace substrings
- Save into generated folders
- Fluent interface
- Generated-file tracking

### `sweepable`

- Variables (`.set`)
- Variable interpolation (`$var`)
- Nested loops (`.foreach`)
- File templating
- Batch generation
- Collection export



## Example

Template file:

```txt
solver cg
dt 0.1
```

Command file:

```txt
.foreach solver cg gmres
    .foreach dt 0.1 0.5
        .read input.txt
        .findLineStartingWith solver
        .replaceBy "solver $solver"
        .findLineStartingWith dt
        .replaceBy "dt $dt"
        .createFolder "output/${solver}_dt_$dt"
        .saveInFolder
    .endforeach
.endforeach
```

Generated files in folders:

```txt
output/cg_dt_0.1/input.txt
output/cg_dt_0.5/input.txt
output/gmres_dt_0.1/input.txt
output/gmres_dt_0.5/input.txt
```



## Command-line language



### General syntax

Commands start with a dot:

```txt
.read input.txt
```

Comments begin with `#`:

```txt
# this is a comment
```

Quoted strings preserve spaces:

```txt
.appendLine "hello world"
```


### Variables

Variables are defined using:

```txt
.set name value
```

Example:

```txt
.set solver gmres
.set dt 0.5
```

Use variables with `$var` (or `${var}`):

```txt
.replaceBy "solver $solver"
```


### Variable Interpolation

Variables are expanded inside arguments.

Example:

```txt
.createFolder "output/$solver"
```

You may also use `${solver}` to avoid ambiguity.

Example:

```txt
.createFolder "output/${solver}_dt_$dt"
```


### Loops

Syntax:

```txt
.foreach variable value1 value2 ...
    ...
.endforeach
```

Example:

```txt
.foreach dt 0.1 0.5 1.0

    .read input.txt
    .findLineStartingWith dt
    .replaceBy "dt $dt"
    .createFolder "output/dt_$dt"
    .saveInFolder

.endforeach
```


### Nested Loops

Nested loops are fully supported.

Example:

```txt
.foreach solver cg gmres

    .foreach dt 0.1 0.5

        .read input.txt

        .findLineStartingWith solver
        .replaceBy "solver $solver"

        .findLineStartingWith dt
        .replaceBy "dt $dt"

        .createFolder "output/${solver}_dt_$dt"
        .saveInFolder

    .endforeach

.endforeach
```



## Commands processor


### `.read`

Load a text file into memory.

```txt
.read input.txt
```


### `.setOutputFilename`

Set output filename.

Default:

```txt
input.txt
```

Example:

```txt
.setOutputFilename config.txt
```



### `.findLineStartingWith`

Select the first line beginning with a prefix.

```txt
.findLineStartingWith dt
```


### `.replaceBy`

Replace the currently selected line.

```txt
.replaceBy "dt 0.5"
```


### `.replaceAllStartingWith`

Replace all matching lines.

```txt
.replaceAllStartingWith solver "solver gmres"
```


### `.replaceInLine`

Replace substring inside selected line.

```txt
.replaceInLine old new
```


### `.insertAfter`

Insert line after matching prefix.

```txt
.insertAfter solver "preconditioner ilu"
```


### `.appendLine`

Append line at end of file.

```txt
.appendLine "# generated automatically"
```


### `.createFolder`

Create output directory.

```txt
.createFolder "output/run1"
```


### `.saveInFolder`

Save current file into current folder.

```txt
.saveInFolder
```

Or directly:

```txt
.saveInFolder "output/run1"
```


### `.silentIfNotFound`

Enable or disable missing-match errors.

Enable:

```txt
.silentIfNotFound
```

Disable:

```txt
.silentIfNotFound false
```


### `.clearCollection`

Clear generated-file collection.

```txt
.clearCollection
```


### `.saveCollection`

Save generated-file collection to text file.

```txt
.saveCollection generated.txt
```

Example output:

```txt
output/cg/input.txt
output/gmres/input.txt
output/dt_0.1/input.txt
```


## Generated File Collection

SFManip automatically tracks every file produced by `.saveInFolder`

The list is stored globally and may be exported using:

```txt
.saveCollection generated.txt
```

This is useful for:

* postprocessing
* batch launching
* archiving
* workflow automation



# C++ API
---

## Basic Usage

```cpp
#include "simuFileManips.hpp"

int main() {

    SFManip sf;

    sf.read("input.txt")
      .findLineStartingWith("dt")
      .replaceBy("dt 0.5")
      .createFolder("output/run1")
      .saveInFolder();

    return 0;
}
```


## Fluent Interface

All methods return:

```cpp
SFManip&
```

allowing chaining:

```cpp
sf.read(...)
  .replaceBy(...)
  .saveInFolder(...);
```


## Main API Methods


### `read`

```cpp
SFManip& read(const std::string& filename)
```

Load file into memory.


### `setOutputFilename`

```cpp
SFManip& setOutputFilename(const std::string& name)
```

Set output filename.


### `findLineStartingWith`

```cpp
SFManip& findLineStartingWith(
    const std::string& prefix
)
```

Select first matching line.


### `replaceBy`

```cpp
SFManip& replaceBy(const char* fmt, ...)
```

Replace selected line using formatted string.

Example:

```cpp
sf.replaceBy("dt %.3f", dt);
```


### `replaceAllStartingWith`

```cpp
SFManip& replaceAllStartingWith(
    const std::string& prefix,
    const char* fmt,
    ...
)
```

Replace all matching lines.


### `replaceInLine`

```cpp
SFManip& replaceInLine(
    const std::string& target,
    const std::string& replacement
)
```

Replace substring inside selected line.


### `insertAfter`

```cpp
SFManip& insertAfter(
    const std::string& prefix,
    const std::string& newLine
)
```

Insert line after matching line.


### `appendLine`

```cpp
SFManip& appendLine(
    const std::string& line
)
```

Append line at end of file.


### `createFolder`

```cpp
SFManip& createFolder(
    const char* fmt,
    ...
)
```

Create and store output folder.


### `saveInFolder`

Stored-folder version:

```cpp
SFManip& saveInFolder()
```

Direct-save version:

```cpp
SFManip& saveInFolder(
    const char* fmt,
    ...
)
```


### `silentIfNotFound`

```cpp
SFManip& silentIfNotFound(bool value=true)
```

Ignore missing matches.



## Static Collection API


### `SFManip::Collection`

Global list of generated file paths.

```cpp
inline static std::vector<std::string> Collection;
```


### `clearCollection`

```cpp
static void clearCollection()
```

Clear generated-file list.



## `saveCollection`

```cpp
static void saveCollection(const char* path)
```

Save collection to text file.

---

# Example: Full C++ Sweep

```cpp
#include "simuFileManips.hpp"

int main() {

    SFManip::clearCollection();

    std::vector<std::string> solvers = {
        "cg",
        "gmres"
    };

    std::vector<double> dts = {
        0.1,
        0.5
    };

    for (const auto& solver : solvers) {

        for (double dt : dts) {

            SFManip sf;

            sf.read("input.txt")
              .findLineStartingWith("solver")
              .replaceBy("solver %s", solver.c_str())

              .findLineStartingWith("dt")
              .replaceBy("dt %.3f", dt)

              .createFolder(
                  "output/%s_dt_%.3f",
                  solver.c_str(),
                  dt
              )

              .saveInFolder();
        }
    }

    SFManip::saveCollection("generated.txt");

    return 0;
}
```


