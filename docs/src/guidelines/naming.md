# General considerations

- When possible, use concise mathematical formulations for equations, making use of convenient ASCII character, like Greek letters, etc.

# Variables

| Long name                             | Variable name                         |
| ------------------------------------- | ------------------------------------- |
| 1st planar coordinate (east-west)     | `x`                                   |
| 2nd planar coordinate                 | `y`                                   |
| Elevation                             | `z`                                   |
| Depth-averaged ice velocity in x      | `vx`                                  |
| Depth-averaged ice velocity in y      | `vy`                                  |
| Depth-averaged ice velocity in z      | `vz`                                  |
| 3D ice velocity in x                  | `vx3D`                                |
| 3D ice velocity in y                  | `vy3D`                                |
| 3D ice velocity in z                  | `vz3D`                                |
| Ice thickness                         | `H`                                   |
| Elevation of ice base                 | `z_base`                              |
| Elevation of ice surface              | `z_surface`                           |

# Functions

Functions should be performing the operations in place. For example `delx!(dudx)`
updates the derivate dudx. We omit the use of a prefix (e.g. `update_delx!`) since
the exclamation at the end of the function name is implicitly understood as an update
of the first argument taken by the function.