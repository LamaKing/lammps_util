# Lammps Util

Collection of scripts to deal with LAMMPS, coming from VASP, CASM and homemade simulations codes.

I'm building this as I go, so take with a pinch of salt.

## Extract Thermo
Parse a LAMMPS output file to extract thermo printout, that usually contain all the interesting stuff.
Define a Python function (to integrate in other scripts or notebooks) and a CLI. Based on Pandas CVS utility and two string flags to decide where to start parsing and where to stop.

Python function returns a list of Pandas Dataframes, one for each run.
CLI either prints to screen all the dataframe or each run to a different file, see help.

