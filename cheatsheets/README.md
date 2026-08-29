# Rockable cheatsheets

One-page reference sheets for `Rockable` and its companion tools, written in
[Typst](https://typst.app) with the
[`boxed-sheet`](https://typst.app/universe/package/boxed-sheet) template.
Each `.typ` compiles to a single landscape A4 page.

| Sheet | Covers |
|---|---|
| `configuration_rockable_cheatsheet` | Structure of a conf-file: timing, neighbor list, strategies, shapes, particles, interactions, interfaces, periodic cell, soft particles, and a minimal input file |
| `physics_rockable_cheatsheet` | Force laws and their parameters, bonds, `Tempo`, global dissipation, body forces, integration schemes, critical time step |
| `companionfiles_rockable_cheatsheet` | `drivingSystem.txt` (controls, servos, periodic loading), `dataExtractors.txt`, pre-processing commands, `postpro`, output files |
| `see_rockable_cheatsheet` | The `see` viewer: CLI, mouse, every key binding, `see.json`, and the other viewers (`seer`, `conftovtk`) |
| `cli_rockable_cheatsheet` | Command lines of every executable, build and CMake options, layout of the repository and of a simulation folder |
| `generator_rockable_cheatsheet` | The `generator` scripting language: output, arithmetic, transformations, shape and packing generation |
| `shapes_rockable_cheatsheet` | Shape-file format, `shapeSurvey`, and the converters `stl2shape`, `tess2shape`, `shape2mesh` |
| `sweepable_rockable_cheatsheet` | `sweepable` and `batchable.sh`: script language, C++ API, and a full parameter study |

## Building

```sh
make          # all the PDFs
make watch    # live preview, e.g. make watch SHEET=physics_rockable_cheatsheet
make clean
```

or, for a single sheet:

```sh
typst compile physics_rockable_cheatsheet.typ
```

The `boxed-sheet` package is fetched by Typst on first use.

## Conventions

`common.typ` holds the shared page setup and the helpers used by every sheet:

- `rockable-sheet` — the page template (landscape A4, four columns). A sheet
  that needs more room lowers `font-size` and `line-skip` rather than spilling
  onto a second page.
- `command(...)` — a keyword, an option value, or a command line.
- `file(...)` — a file name.
- `flag(...)`, and the shorthands `periodic-flag`, `soft-flag`,
  `boundary-flag`, `ftcorr-flag` — mark a feature that needs a CMake option to
  be enabled. `rockable -b` prints which ones a binary was built with.
- `note[...]` — a caveat or an aside.
- `kwtable(...)` — a two-column keyword table. Use a bulleted list instead when
  the keys are long, otherwise the description column gets squeezed.

The content is drawn from the source rather than from the documentation, so a
few entries deliberately differ from `sphinxdoc`, in particular the
`tranformation:` keywords of `generator` and the `MinkowskiRadius` keyword of
`tess2shape`.
