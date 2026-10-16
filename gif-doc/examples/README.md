# Animated GIFs of the examples

These scripts make the `animation.gif` files of the folders of `examples`,
shown in the readmes and in the galleries of `README.md` and
`examples/README.md`.

They need Rockable built in `BUILD` (with `see`, `generator` and
`stl2shape`, as done by `install_rockable.sh`), and `ffmpeg`, `gifsicle` and
ImageMagick (`magick`). Python is not used.

## Usage

```sh
cd gif-doc/examples
sh build_see_auto.sh            # once, and after each change of see.cpp
sh all_gifs.sh                  # all the examples (about 45 minutes)
sh all_gifs.sh helloworld       # or only some of them
sh all_gifs.sh install          # copies work/*.gif into the examples
```

A `see` window opens on the screen while the frames are rendered: leave it
alone until it closes by itself.

Everything is written in `work/`: one folder per example (the copy of the
example where it was run, with its conf-files and its frames),
`work/<name>.gif` and a contact sheet `work/<name>.png` (four frames side by
side, to check the framing and the colours before looking at the GIF).

## The scripts

| File                | Role                                                                                         |
|---------------------|----------------------------------------------------------------------------------------------|
| `autoshot.inc`      | A function added to `see`: renders the conf-files one after the other into `frameNNNN.tga`, then exits |
| `build_see_auto.sh` | Inserts `autoshot.inc` into a copy of `src/Apps/see.cpp` and compiles `see_auto` with the flags of `BUILD` |
| `run_example.sh`    | Copies an example into a work folder, sets `interConf`/`tmax` (shortened runs), runs `rockable` |
| `render.sh`         | Runs `see_auto` with a `see.json` (camera, colours), then `ffmpeg` (palette) and `gifsicle` (compression) |
| `sheet.sh`          | Contact sheet of four frames                                                                 |
| `autocam.sh`        | Writes a `see.json` whose camera frames the particles of a conf-file (a starting point for a new example) |
| `all_gifs.sh`       | The recipe of each example: camera, shortened run, colour range, size of the GIF              |

`see_auto` is driven by environment variables (set by `render.sh`):
`AUTOSHOT_STEP` (conf-files between two frames), `AUTOSHOT_LAST` (last
conf-file), `AUTOSHOT_CMIN` and `AUTOSHOT_CMAX` (colour range fixed for the
whole animation, otherwise `see` rescales it on each conf-file).

## Adding an example

1. Run it with `run_example.sh` (shorten it, and choose `interConf` so as to
   get 50 to 100 conf-files).
2. Find a camera: `sh autocam.sh work/<name>/conf0 2.0 0.5 0.4 1 > work/<name>.json`
   (distance factor, then viewing direction), or open `see` in the work
   folder, move the view and save it in `see.json` (key `s`).
3. `sh render.sh work/<name> work/<name>.json work/<name>.gif` and
   `sh sheet.sh work/<name> work/<name>.png`; adjust.
4. Write the recipe in `all_gifs.sh` (and in `dest` if the folder is nested).

Keep the GIFs light (below about 300 KB): reduce `WIDTH`, `COLORS`, the number
of frames (`LAST`, `STEP`, `FRAMES`), or increase `LOSSY`. Scenes with many
textured particles (polyhedra, meshes) are the heaviest; `DITHER=none` helps.
