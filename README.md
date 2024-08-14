# Tug of War
This is some code to graph rupture probabilities of 2 bonds in a _Tug of War_ simulation.

### `Main.m`
This calls functions from the other scripts and plots the corresponding graphs.

### `Potentials.m`
A set of utility functions provided for `Main.m`, including bond rupture energies.

### `antigen-extract-dynamics.m`
Performs a random walk simulation on the energy landscape, using _Euler Muryama (Stochastic)_ time stepping.

<img src="https://github.com/johoesel/Tug-of-war/blob/main/Extraction_vs_Diffusion.png" width="600">
