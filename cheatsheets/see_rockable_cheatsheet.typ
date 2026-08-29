#import "common.typ": *

#show: rockable-sheet.with(
  title: [`see`, the 3D visualiser for `Rockable`],
  write-title: true,
)

= Command line

#concept-block(body: [
  ```sh
  see                # opens conf0
  see conf10         # opens the configuration number 10
  see input.txt      # opens a conf-file by name
  see -t traj.txt    # loads a trajectory file
  see -v 4           # verbose level of the Rockable core
  ```

  `see` opens the conf-file in *interactive mode*: no data extractor is run and
  nothing is written back, so the configurations are safe to browse.
])

= Files picked up automatically

#concept-block(body: [
  #kwtable(
    file("see.json"), [Every display parameter and the camera. Loaded at
      start-up if present; written with #command("J"), reloaded with
      #command("j").],
    file("<shapeFile>.rmsh"), [Skin meshes next to the shape library (produced
      by `shape2mesh --rmsh`). When present, the particles are drawn with their
      meshes instead of the overlapping primitives.],
    file("traj.txt"), [Trajectories to overlay. The name is settable with `-t`.],
    file("probe.txt"), [A measurement box: `min` (3 reals), `max` (3 reals) and
      the number of Monte-Carlo steps. Toggled with #command("!"), evaluated
      with #command("@").],
    file("drivingSystem.txt"), [Read like in a computation, so the driven bodies
      are shown consistently.],
  )
])

= Mouse

#concept-block(body: [
  #kwtable(
    command("left + move"), [Rotate.],
    command("shift + left + move"), [Pan.],
    command("ctrl + left + move"), [Zoom (the middle button also works).],
    command("shift + middle"), [Select a particle (same as #command("*")).],
  )
  The key #command("h") shows this list in the window.
])

= Navigation and interface

#concept-block(body: [
  - #command("+ / -") ~Go to the next / previous conf-file.
  - #command("g") ~Open another conf-file; the number is asked *in the terminal*
    (a negative number escapes).
  - #command("=") ~Fit the view to the whole scene.
  - #command("w") ~Set the vertical of the view along gravity.
  - #command("*") ~Select the particle under the mouse; #command("ESC")
    deselects.
  - #command("i") ~Print an info panel: the whole system (number of bodies,
    `t`, `dt`, `tmax`, the Verlet distances and the ratio $Delta t_c \/ Delta t$),
    or the selected particle (group, cluster, position, velocities, quaternion,
    shape, mass and inertia).
  - #command("p") ~Edit the selected body, in the terminal.
  - #command("x") ~Print the bounding limits of the scene, at full precision in
    the terminal.
  - #command("h") ~Show the mouse usage. #command("k") ~Show the key bindings.
  - #command("space") ~Clear the text at the bottom of the window.
  - #command("arrow up / down") ~Change the number of text lines at the bottom.
  - #command("q") ~Quit.
])

= What is displayed

#concept-block(body: [
  - #command("n") ~Bodies on/off.
  - #command("d") ~Driven bodies (the solid boundaries) on/off.
  - #command("a / A") ~Decrease / increase the transparency of the driven
    bodies (step of 0.05).
  - #command("e / E") ~Decrease / increase the transparency of the free bodies.
  - #command("o") ~Oriented Bounding Boxes on/off.
  - #command("O") ~Enlarge the OBBs by half of the secure distance.
  - #command("f") ~Force vectors at the contacts.
  - #command("F") ~Centre-to-centre normal forces.
  - #command("l") ~Local frames at the contacts.
  - #command("m") ~Contacts coloured by sub-interaction type.
  - #command("P") ~Periodic cell #periodic-flag.
  - #command("!") ~Probe box on/off (needs #file("probe.txt")).
  - #command("@") ~Compute the solid fraction inside the probe, by Monte-Carlo
    sampling; the result is printed in the terminal.
  - #command("b") ~Background colours on/off.
  - #command("y") ~Faster display: the shapes are drawn without their
    Minkowski thickness. Handy on large samples.

  #note[#command("v") toggles the `show_velocities` flag, but the drawing of the
  velocity vectors is currently commented out in the code.]
])

= Colours

#concept-block(body: [
  - #command("0") ~All the particles in the same colour.
  - #command("1") ~Cyclic colours, one per particle.
  - #command("2") ~Coloured by velocity magnitude.
  - #command("3") ~Coloured by shape name.
  - #command("c") ~Colour bar on/off.
  - #command("r") ~Toggle the automatic rescaling of the colour range. When
    switched off, the range is taken from `colorRangeMin` and `colorRangeMax` of
    #file("see.json"); the current values are echoed in the terminal.
])

= Screenshots

#concept-block(body: [
  - #command("z") ~One screenshot, into `oneshot.png` (or `oneshot.tga` when
    `see` was built without `libpng`).
  - #command("Z") ~One screenshot per conf-file, from the current one to the
    last readable one: `shot<N>.png`. #note[There is no way to interrupt this
    loop, so check the number of conf-files first.]

  A quick animation from the series:
  ```sh
  ffmpeg -framerate 25 -i shot%d.png -pix_fmt yuv420p movie.mp4
  ```
])

= Debugging aid

#concept-block(body: [
  - #command("%") ~Rebuild the neighbor list, then run 5000 velocity-Verlet
    steps on the configuration being displayed. Nothing is saved; this is only
    meant to see where a configuration is heading.
])

= `see.json`

#concept-block(body: [
  Written by #command("J"), re-read by #command("j") or at start-up. Handy to
  reproduce exactly the same view over several samples, or to set a colour range
  that the interface cannot reach.

  ```json
  {
    "ParticleColor": [207, 174, 85],
    "alpha_fixparticles": 0.15,
    "alpha_particles": 1.0,
    "camera": {
      "center": [0.0, 0.0, 0.0],
      "eye": [-0.715, 0.400, 0.573],
      "up": [0.313, 0.916, -0.250],
      "view_angle": 45.0,
      "znear": 0.1,
      "zfar": 2.194
    },
    "colorMode": 0,
    "colorRangeMin": 0.0,
    "colorRangeMax": 1.0,
    "rescaleColorRange": 1,
    "enlarged_obb": 0,
    "show_background": 0,
    "show_colorBar": 1,
    "show_driven": 1,
    "show_forces": 0,
    "show_interFrames": 0,
    "show_interTypes": 0,
    "show_keybinds": 0,
    "show_normal_forces": 0,
    "show_obb": 0,
    "show_particles": 1,
    "show_periodicBox": 1,
    "show_probe": 0,
    "show_traj": 0,
    "show_velocities": 0,
    "window": { "width": 768, "height": 768 }
  }
  ```

  #note[`show_traj` has no key binding: switch it on here, with a
  #file("traj.txt") next to the conf-files.]
])

= The other viewers

#concept-block(body: [
  #inline("seer --- Dear My Seer")
  The successor of `see`: same input format, but built on *Dear ImGui* and
  *SDL* instead of `freeglut`/X11, which makes it much easier to maintain on
  macOS. Most of what `see` puts behind a key is a widget of a control panel,
  and it adds slicing planes and a richer information panel.

  - Selection: #command("ctrl + left click") on a particle; #command("ctrl + left click")
    on the background deselects.
  - Built when `ROCKABLE_COMPILE_SEER` is on. Its sibling *Dear My Shape*
    (`src/Apps/DearMyShape`, binary `shape`) does the same for shape libraries,
    but has its own `Makefile`.

  #note[`seer` reads the Rockable conf-file format, so it also displays the
  output of other codes that write it --- `ExaDEM` for instance.]

  #inline("conftovtk --- ParaView")
  ```sh
  conftovtk
  ```
  Converts *every* `conf*` file found in the current folder into VTK legacy
  files, ready for ParaView:
  #kwtable(
    file("Shepres<N>.vtk"), [The spherical particles.],
    file("Polyr<N>.vtk"), [The polyhedral particles.],
    file("PolyrSpheres<N>.vtk"), [Their vertex spheres.],
    file("Forces<N>.vtk"), [The contact forces.],
    file("OBB<N>.vtk"), [The oriented bounding boxes.],
  )
  Use it for the figures you want to compose precisely, or for the fields `see`
  cannot draw.
])
