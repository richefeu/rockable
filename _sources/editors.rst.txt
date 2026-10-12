.. _editors:

Editing the input files: ``rockedit`` and ``confedit``
======================================================

Two editors know the language of the Rockable input files (conf-files,
shape files, ``drivingSystem.txt``, the command files of ``postpro``): they
colour the keywords, **show the documentation of a keyword**, and insert
ready-made blocks (snippets). They are the quickest way to write an input file
without keeping this documentation open next to it.

``rockedit``
------------

A small editor that runs in the terminal, with no dependency, so that it can
also be used over ``ssh`` on a cluster. It is built and installed with
Rockable.

.. code-block:: sh

   rockedit input.txt

The keys are those of ``nano``, and the bottom bar always shows them:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Key
     - Action
   * - ``Ctrl+D``
     - Show (or hide) the documentation of the word under the cursor, or of
       the selection. The pane follows the cursor: move onto ``tmax`` and it
       explains ``tmax``.
   * - ``Ctrl+P``
     - Insert a snippet (a soft particle, a periodic cell, a spring joint...);
       type to filter the list.
   * - ``Ctrl+O`` or ``Ctrl+S``
     - Save.
   * - ``Ctrl+X``
     - Quit (offers to save a modified file).
   * - ``Ctrl+W`` / ``Ctrl+N``
     - Find / find next.
   * - ``Ctrl+G``
     - The full list of keys.

The file ``prepro/rockedit/README.md`` describes all the keys.

``confedit``
------------

The graphical counterpart of ``rockedit`` (it needs FLTK), with the same
colouring, documentation and snippets.

.. code-block:: sh

   confedit input.txt

Where the documentation comes from
----------------------------------

Both editors read the file ``prepro/common/rockable.lang``, which describes
every keyword of the language, its arguments and its meaning, and holds the
snippets. A copy is placed next to the executables when they are built.
Documenting a new keyword in this file is enough for both editors to show it,
without recompiling them.
