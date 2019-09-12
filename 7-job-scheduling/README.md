# Chapter 7. Job scheduling

The simulation results of Section 7.3 can be obtained and plotted by running the following two files:
- ``compute.sh``: Compute the simulation results (takes a few hours on a commodity computer). Compile the C source files contained in the ``sources/`` folder and calls the obtained executables. The results are stored in the folder ``data/``.
- ``plot.ipynb``: Plot the simulation results and compare them with the exact results predicted by the model of Section 8.2.

The folder ``sources/`` gathers the source files for the simulations, written in C language:
- ``toy-*.c``: Simulation files for the toy example of Section 7.3.1.
- ``random-*.c``: Simulation files for the random job assignment of Section 7.3.2.
- ``utils.h`` and ``utils.c``: Libraries, constants, and random sampling functions used in the C programs.
