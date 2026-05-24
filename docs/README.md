# Overview

This explains how to generate the documentation for ImpactX, and contribute to it.

## Generating the documentation

### Installing the requirements

`cd` into the ImpactX `docs/` directory.
Install the Python requirements for compiling the documentation with **one** of the following package managers:

- pip: `pip install -U -r requirements.txt`
- uv: `uv pip install -U -r requirements.txt`
- conda: `conda create -n readthedocs -f conda.yml && conda activate readthedocs`
- spack: `spack env activate -d . && spack install`

### Compiling the documentation

`cd` into this directory and type
```
make html
```
You can then open the file `build/html/index.html` with a standard web browser (e.g. Firefox), in order to visualize the results on your local computer.

### Cleaning the documentation

In order to remove all of the generated files, use:
```
make clean
```
