# hydro v0.3
#### Last Update 10/28/2026
This repository holds the hydrodynamics code that I have been developing as a graduate student at Clemson University (Fall 2025-Present). This code currently has the following implemented:
- 1D Hydrocode - Cartesian, Spherical (In progress!)
- 1D MHD - Cartesian
- 2D Hydrocode - Cartesian

### Code Structure
#### main.rs
This is where the simulation is initialized and run. The performance metrics are output by this file while the simulation is run. This file is responsible for writing txt files after each checkpoint for plotting purposes.

#### math_func.rs
This file holds many of the formulas required to do the calculations necessary for this simulation. This includes physical formulas such as the formula for the sound speed of the gas. The formulas to compute the timestep and to perform piecewise reconstruction are here too.

#### sim_func.rs
In this file, the functions that define the grid and advance the timestep are defined. This file also has the function that writes the checkpoint files. At the top of this file the data structures used throughout the entire code.

#### test_conds.rs


#### math_functions.py


#### plotting_functions.py


#### simple_animator.py


### Planned Features