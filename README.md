[![Anaconda-Server Badge](https://anaconda.org/dbankmann/pygelda/badges/version.svg)](https://anaconda.org/dbankmann/pygelda) [![Anaconda-Server Badge](https://anaconda.org/dbankmann/pygelda/badges/latest_release_date.svg)](https://anaconda.org/dbankmann/pygelda)
# Summary
pygelda is a python interface for the FORTRAN solver GELDA. The interface is currently in alpha state and only supports strangeness-free problems.

See `help(pygelda.Gelda)` for rudimentary help on the interface.

pygelda ships with the current version of GELDA v1.1.1, for which the license in `GELDA/LICENSE` applies. It is available from https://doi.org/10.1137/S1064827595286347.

# Installation

Linux users can use the provided conda package `pygelda` which can be installed by
```
conda install -c dbankmann pygelda
```
