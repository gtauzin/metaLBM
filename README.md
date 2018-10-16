## Description

metaLBM is a C++11 header-only template Computational Fluid Dynamic software
based on the Lattice Boltzmann Method (LBM). metaLBM is meant to support multiple lattices
on highly parallel architecture. While it can be used to perform various kind
of flow simulations, it has been developped with the aim of exploring turbulent flow
simulation and turbulence modelling within LBM.The library is hybrid parallel with
MPI for distributed memory and OpenMP threads for shared memory. It is under active
development as of April 2017 and is currently being ported on multi-GPU architectures.

The [Lattice Boltzmann Method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods)
is a meso-scale approach to the simulation of fluid dynamics. Instead of solving the
Navier–Stokes equations, the discrete Boltzmann equation is solved to simulate the flow
of a Newtonian fluid with collision models such as Bhatnagar–Gross–Krook (BGK). By
simulating streaming and collision processes across a limited number of particles, the
intrinsic particle interactions evince a microcosm of viscous flow behavior applicable
across the greater mass.

The project provides a general framework for experimentation with simulations
of turbulent flows. It has an easy mechanism allowing to use various forcing schemes,
forces, and boundary conditions. The code currently runs on several medium size clusters
Stromboli at Bergische Universitaet Wuppertal (BUW) and Newturb at Università degli Studi
di Roma "Tor Vergata" and on larger supercomputers such as Galileo at CINECA in Bologna.
It is expected to run soon on JURECA and JURON at the Juelich Supercomputing Center (JSC).

## Overview

As of April 2017, metaLBM support the following features:

**Supported lattices**
- D1Q3
- D2Q5
- D2Q9
- D3Q15
- D3Q19
- D3Q27
- ... potentially any rectangular-shaped lattice

**Supported LBM algorithm**
- Fused collide and stream pull
- Fused collide and stream push (upcoming for selected lattices)

**Single-node performance optimization**
- SoA and AoS memory layout
- Automatic OpenMP vectorization (upcoming)
- OpenMP multi-threading

**Multi-node performance optimization**
- 1D domain partionning
- 2D and 3D domain partitionning (upcoming)

**Supported collision operators**
- SRT BGK
- SRT Entropic
- MRT Entropic (upcoming)

**Supported forcing schemes**
- Shan-Chen
- Guo
- Exact-Difference Method

**Supported forces**
- Constant
- Sinusoidal
- Kolmogorov
- Turbulent forcing on a shell of wavenumbers with time-independent phase (upcoming)
- Turbulent forcing on a shell of wavenumbers with time-dependent phase (upcoming)

**Supported boundary conditions**
- Periodic BC
- Half-way bounceback (upcoming)
- Entropic bounceback (upcoming)

**Supported I/O**
- serial VTK (ascii, binary upcoming)
- HDF5 (upcoming)

## Project Organisation

The project structure follows a common C++ layout within the [solver](solver) directory.

- The [test](test) directory contains source code for unit tests.
- The [include](include) directory contains headers-only library dependencies.
- The [src](src) directory contains all source files.
- The [log](log) directory contains log files generated when ENABLE_LOG
is defined.
- The [doc](doc) directory contains documentation generating code
with doxygen.

The [src](src) directory contains several subdirectories which are aiming
at exploiting different type of parallelism. They are included if the corresponding option
is passed to CMake.


Output files are dumped in the [output](output) directory.
- The [output/outputVTR](output/outputVTR) directory contains all generated VTK
files.

## Building metaLBM

Provided a working installation of CMake(>= 3.8) is present, `cd` to project
root and setup build directory

```shell
mkdir build
cd build
cmake ..
```

If the system's compiler is too old and does not support `C++11`, you can specify
a custom compiler:

```shell
cmake .. -DCMAKE_CXX_COMPILER=/path/to/compiler
```

metaLBM requires a `MPI 3.0` compliant library implementation is needed (eg
`OpenMPI`, `MPICH`). If the library is already installed, `CMake ` usually
finds the library without problems. In case multiple `MPI` libraries are installed, do
To build tests, you will need a working intallation of
[GoogleTest](https://github.com/google/googletest). You may have to set the
environment variable `GTEST_ROOT` to the location of the `GTest`'s root directory
before running `cmake`:

```shell
export GTEST_ROOT="path/to/gtest/root"
```

To build all tests do:

```shell
make tests
```

If you want to add unit tests, just add a source file in [test/](test/),
edit the correspondng `CMakeLists.txt` files to include your new targets and rebuild the
project `cmake ..` from your [build](build) directory. There are
predeclared functions which make it easy to add new targets.

metaLBM supports forcing in Fourier space, which is based on [FFTW](http://www.fftw.org)
library and therefore you will need a working installation. You may have to set the
environment variables `FFTW_INCLUDE_DIR` to the location of the `FFTW`'s include directory,
and `FFTW_LIBRARY_DIR` to the location of the `FFTW`'s lib directory before running `cmake`:

```shell
export FFTW_INCLUDE_DIR="path/to/fftw/include"
export FFTW_LIBRARY_DIR="path/to/fftw/lib"
```

To generate documentation for the project, you need a working `doxygen`
installation. Run from your [build](build) directory:

```shell
make doc
```

Documentation will be generated in a subfolder [doc](doc) and the webpage
can be found in [doc/html/index.html](doc/html/index.html).

Check the [Wiki](https://gitlab.com/rooknrowl/metaLBM/wikis/home) for
additional tips on how to run a simulation.


## License
Please see [LICENSE](LICENSE).


## Contributing
Please see [CONTRIBUTING.md](CONTRIBUTING.md)
