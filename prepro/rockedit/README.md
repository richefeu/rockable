rockedit
========

A small terminal editor for the Rockable input files (conf, `shapes.txt`,
`drivingSystem.txt`, postpro inputs). It is the terminal counterpart of
`prepro/confedit`: same syntax colouring, same inline documentation, but light
enough to be used over ssh on a cluster.

No dependency: termios and ANSI escape sequences only, so a plain C++17 compiler
is all it takes.

Build
-----

It is built with the rest of Rockable (the top-level `CMakeLists.txt` carries the
compiler, the standard and the install prefix; `prepro/CMakeLists.txt` pulls this
directory in).

Standalone, when all you have is a compiler on a cluster:

~~~bash
make
./rockedit ../../examples/DominoFun/input.txt
~~~

Keys
----

Nano-like; the bottom bar always shows them, and `^G` opens the full list.

| Key | Action |
|---|---|
| `^O` (or `^S`) | save; asks for a name when there is none |
| `^X` | quit (offers to save when the buffer is modified) |
| `^W` / `^N` | find / find next, case-insensitive, wraps around |
| `^L` | go to a line number |
| `^D` | toggle the documentation of the selection, or of the word under the cursor |
| `^P` | insert a snippet (type to filter the list) |
| `Shift` + arrows | select; `Shift` with `Home`, `End`, `PgUp`, `PgDn` too |
| `^K` / `^U` | cut the selection, or the whole line when there is none / paste |
| `^Z` / `^Y` | undo / redo |
| `^A` / `^E` | start / end of the line (`Home` and `End` too) |
| `Ctrl+Up` / `Ctrl+Down` | start / end of the file |

The documentation pane opened by `^D` follows the cursor: move onto `tmax` and it
explains `tmax`. The status bar hints `^D documents '<word>'` whenever the word
under the cursor is documented and the pane is closed.

Selection
---------

Hold `Shift` and move: the span between where you started and the cursor is shown
in reverse video, and the status bar counts it. A movement without `Shift` drops
the selection.

Typing, `Return`, `Tab`, `Backspace`, `Delete` and `^U` all replace what is
selected. `^D` documents the selected text rather than the word under the cursor,
the way `confedit` does on its own selection — so selecting `gravity` and pressing
`^D` explains `gravity`.

`^K` cuts the selection and `^U` pastes it back, at the cursor. With no selection
`^K` still cuts the whole current line, and `^U` then puts it back as a whole
line rather than in the middle of another one. There is no separate copy key:
cut with `^K`, press `^U` once to put it back where it was, then `^U` again
wherever you want the copy.

`Shift+PgUp` and `Shift+PgDn` may be swallowed by the terminal itself, which
often binds them to its own scrollback.

Where the colours and the documentation come from
-------------------------------------------------

Everything the editor knows about the Rockable language sits in
`prepro/common/rockable.lang`: 117 keywords, 55 types, 116 documentation entries and 23
snippets. `confedit` reads the very same file, so the two editors cannot drift
apart, and documenting a new keyword needs no recompilation.

The file is looked up, in this order:

1. `$ROCKABLE_LANG`
2. the current directory
3. next to the executable, then in its three parent directories (CMake drops a
   copy next to the binaries at build time)
4. `~/.rockable`

Each of those is tried both directly and through a `common/` subdirectory, which
is where the file sits in the source tree.

Its format is described in its own header. In short:

~~~
[keyword] tmax
doc: tmax [(double) value]

  Maximum time (at the end of a simulation)

[type] velocityVerlet

[snippet] Periodic cell
usePeriodicCell 1
h _xx_ _xy_ _xz_ _yx_ _yy_ _yz_ _zx_ _zy_ _zz_
~~~

A `[keyword]` is coloured blue, a `[type]` cyan, and a token starting with `/`,
`#` or `!` turns the rest of the line green — exactly the rule `confedit` uses.

Notes
-----

* The whole file is held in memory as a vector of lines and only the visible
  ones are rendered, so opening the 11 MB `examples/BigForTests/input.txt` is
  instantaneous.
* Bytes are treated as columns: the Rockable input files are ASCII.
* `^C` does not quit (signals are disabled in raw mode); use `^X`.

Source layout
-------------

| File | Role |
|---|---|
| `src/main.cpp` | the editor: layout, keys, colouring |
| `src/text_buffer.hpp` | lines, cursor, editing, undo/redo, search |
| `src/terminal.hpp` | raw mode, key decoding, frame output |
| `../common/rockable_lang.hpp` | reads `rockable.lang`; shared with `confedit` |
