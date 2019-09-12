# Chapter 8. Load balancing

The simulation results of Section 8.3 can be obtained and plotted by running the following files:
- ``compute.sh``: Compute all results concerning the dynamic load-balancing algorithm, except for those of the paragraph "Number of tokens" in Section 8.3.1 (takes a few hours on a commodity computer). Compile the C source files contained in the ``sources/`` folder and calls the obtained executables.
- ``dynamic-onetype-exact.ipynb``: Compute the performance of the dynamic load-balancing algorithm when the number of tokens increases, as described in the paragraph "Number of tokens" of Section 8.3.1.
- ``static-exact.ipynb``: Compute the performance of the static load-balancing policies in both scenarios.
- ``plot.ipynb``: Plot the results in a Jupyter Notebook.

All results are stored in the folder ``data/``.

The folder ``sources/`` gathers the source files for the exact and simulation results, written in C language, under the dynamic load-balancing algorithm:
- ``dynamic-exact.c``: Exact results predicted by the model of Section 8.2.
- ``dynamic-exp.c``: Simulation results under exponentially-distributed job sizes.
- ``dynamic-hyperexp.c``: Simulation results under hyperexponentially-distributed job sizes.
- ``utils.h`` and ``utils.c``: Libraries, constants, and random sampling functions used in the C programs.

All obtained data files are stored in the folder ``data/``.
