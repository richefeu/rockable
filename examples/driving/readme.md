## 2. Driving a body

### Goal

A plate rotates at one turn per second, like a flipper, and hits three grains
thrown towards it. The rotation of the plate is imposed in a separate file,
`drivingSystem.txt`: this is how a boundary is moved in Rockable.

### Notions

- **Driven bodies**: the first `nDriven` particles of the list.
- The file **`drivingSystem.txt`** and its `Control` entries, which impose a
  velocity, a force or a moment on one component of a driven body
  ([drivingSystem](https://richefeu.github.io/rockable/drivingSystem.html)).
- The **servo-controllers** (`Servo`), which recompute the controls at every
  step ([servoControllers](https://richefeu.github.io/rockable/servoControllers.html)).
- The force law `Avalanche`, an alternative to `Default`
  ([forceLaws](https://richefeu.github.io/rockable/forceLaws.html)).

### Files

| File                | Content                                                    |
|---------------------|------------------------------------------------------------|
| `input.txt`         | The conf-file: a plate (driven) and three grains            |
| `drivingSystem.txt` | What is imposed on the plate                                |
| `SphereAndPlan.shp` | The shapes `Plan`, `Sphere` and `Rice`                      |

### Run it

```sh
rockable input.txt   # about 1 second
see
```

### What happens

In `input.txt`, `nDriven 1` makes the first particle, the plate, a driven body.
The three grains start 0.5 m above it, with a velocity `(3.69, -3.29, 0)`
m/s, that is about 5 m/s towards the plate. The second elongated grain has a
homothety of `0.5`: it is half the size of the first one.

`drivingSystem.txt` holds a single line:

```
Control _zrot_Vel_ 0 6.28
```

which imposes on body `0` an angular velocity of 6.28 rad/s about the `z` axis.
The plate turns about its centre, and makes one full turn during the second of
simulation. The sphere is hit at about `t = 0.1 s` and thrown upwards at about
15 m/s. The two elongated grains are swept along by the plate, then thrown
away near the end.

The general form of a control is

```
Control <type> <body> <value>
```

| `<type>`                                   | Imposed                                    |
|--------------------------------------------|--------------------------------------------|
| `_x_Vel_`, `_y_Vel_`, `_z_Vel_`            | a velocity component                       |
| `_xrot_Vel_`, `_yrot_Vel_`, `_zrot_Vel_`   | an angular velocity component              |
| `_x_For_`, `_y_For_`, `_z_For_`            | a force component                          |
| `_xrot_Mom_`, `_yrot_Mom_`, `_zrot_Mom_`   | a moment component                         |
| `_xyzrot_Vel_`, `_xyzrot_Mom_`             | the angular velocity, the moment (3 values) |

Three rules are worth remembering:

- A component of a driven body that is not controlled **does not move**. Here,
  the plate keeps its position and only rotates about `z`. Without the file
  `drivingSystem.txt`, the plate would not move at all, which is how fixed
  walls are made.
- A component driven by a force is integrated like a free body, but **gravity
  is not applied** to it: the imposed force is the only external force.
- `drivingSystem.txt` is **read again** each time a conf-file is saved: it can
  be edited while the computation runs.

### What to look at

- In `see`, follow the rotation of the plate with `+`, and the trajectories of
  the grains.
- The angular velocity of the plate is the 3rd component of `vrot` in its line
  of the conf-files; it stays at `6.28`.

### Exercises

1. Reverse the rotation (`-6.28`), or slow it down (`3.14`).
2. Replace the content of `drivingSystem.txt` by `Control _y_Vel_ 0 0.5`: the
   plate no longer rotates, and rises at 0.5 m/s.
3. Put both lines in the file: the plate rises **and** rotates.
4. Replace the content by `Control _y_For_ 0 0`: the plate is now free along
   `y`, with no force imposed and no gravity. It is pushed down by the grains
   that land on it.
5. Make the plate vibrate with a servo: `Servo shaker 0  0 1 0  0.02 5`
   (body 0, direction `(0, 1, 0)`, amplitude 0.02 m, frequency 5 Hz).
6. Delete `drivingSystem.txt`: the plate stays fixed, and the grains bounce on
   it.
7. Start a long computation (set `tmax 5.0`), and edit `drivingSystem.txt`
   while it runs: the change is taken into account at the next conf-file.

### Going further

- [`isoComprSpheres`](../isoComprSpheres): six walls driven by a servo to
  compress a sample.
- All the controls and servos:
  [drivingSystem](https://richefeu.github.io/rockable/drivingSystem.html) and
  [servoControllers](https://richefeu.github.io/rockable/servoControllers.html).
