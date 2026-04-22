.. _Visualisation:

Visualisation tools
===================

Application `see`
-----------------

Normally, the application ``see`` has been built as the same time than ``rockable``. 
The application ``see`` needs ``freeglut``, the simplest way to use openGL and display 3D things.
``see`` allows for the visualisation of conf-file (including an input-file).


**Usage**

Here is the console output when typing ``see --help``:

.. code-block:: text

   USAGE:
   
      ../../BUILD/see  [-v <int>] [-t <string>] [--] [--version] [-h]
                       <conf-file>
   
   
   Where:
   
      -v <int>,  --verbose <int>
        Verbose level
   
      -t <string>,  --traj <string>
        Name of a trajectory file
   
      --,  --ignore_rest
        Ignores the rest of the labeled arguments following this flag.
   
      --version
        Displays version information and exits.
   
      -h,  --help
        Displays usage information and exits.
   
      <conf-file>
        Name of the conf-file
   
   
      Visualisation of Rockable simulations


So, to visualise a conf-file you can use its name

.. code-block:: sh

   see input.txt

or simply invoke ``see`` with no other argument, and the default one will be ``conf0``.



Application `seer` (Dear My Seer)
---------------------------------

.. image:: images/GPT_DearMySeer_logo_white.png
   :width: 250px
   :align: center

**Dear My Seer** is the visualization application that replaces ``see`` in **Rockable**.
The old ``see`` application uses ``freeglut`` and thus **X11**, which made maintenance difficult, especially on Apple computers. The console-name ``seer`` stands for **see-rockable**, and since this new version uses **Dear ImGui**, it has been named **Dear My Seer** as a playful wordplay.

.. note::

   If you were using ``see``, you will only need to slightly change your habits by adding an ``r`` at the end of ``see``.

Thanks to **Dear ImGui**, many features are now accessible through a modern and intuitive interface, using the mouse for a smoother experience.

**Dear My Seer** is based on **Rockable**'s input format, meaning it can be used with other DEM applications.
All you need to do is provide the configuration file (conf-file) in the correct format (this is what **ExaDEM** already does).

