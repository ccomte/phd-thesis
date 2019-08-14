# Chapter 7. Job scheduling

The simulation results of Section 7.3 can be obtained and plotted by running the following two files:
- ``compute.sh``: Compute the results (takes a few hours on a commodity computer). Compile the C source files contained in the ``sources/`` folder and calls the obtained executables.
- ``plot.ipynb``: Plot the results in a Jupyter Notebook.

The ``sources/`` folder gathers the source files for the simulations, written in C language:
- ``toy-*``: Simulation files for the toy example of Section 7.3.1.
- ``random-*``: Simulation files for the random job assignment of Section 7.3.2.
- ``utils.h`` and ``utils.c``: Libraries, constants, and random sampling functions used in the C programs.
