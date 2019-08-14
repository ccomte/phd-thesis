# Chapter 5. Job scheduling

The simulation results can be obtained and plotted with the following two files:
- ``compute.sh``: Run the simulations (takes a few hours on a commodity computer).
- ``plot.ipynb``: Plot the obtained results.

The simulations are written in C language. All source files are in the ``sources/`` folder:
- ``toy-*``: Simulation files for the toy example of Section 7.3.1.
- ``random-*``: Simulation files for the random job assignment of Section 7.3.2.
- ``utils.h`` and ``utils.c``: Libraries, constants, and random sampling functions used in the C programs.
