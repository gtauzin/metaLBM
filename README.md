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
di Roma "Tor Vergata" and on larger supercomputers such as Galileo at CINECA in Bologna. It is expected to run soon on JURECA and JURON at the Juelich Supercomputing Center (JSC).

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

**Supported collision operators**
- BGK
- Entropic BGK

**Supported forcing schemes**
- Shan-Chen
- Guo
- Exact-Difference Method

**Supported I/O**
- VTK (serial)
- HDF5 (upcoming)

## Project Organisation

The project structure follows a common C++ layout within the [solver](solver) directory.

- The [solver/test](solver/test) directory contains source code for unit tests.
- The [solver/include](include) directory contains headers-only library dependencies.
- The [solver/src](solver/src) directory contains all source files.
- The [solver/doc](solver/doc) directory contains documentation generating code
with doxygen.

The [solver/src](solver/src) directory contains several subdirectories which are aiming
at exploiting different type of parallelism. They are included if the corresponding option
is passed to CMake.

- The [solver/src/core](solver/src/core) directory contains core source code describing the framework
and the physics.
- The [solver/src/omp](solver/src/omp) directory contains a serial or shared memory
parallelisation version of the LBM algorithm. Activate it by passing the flags `-D_SERIAL=ON`
or `-D_OMP=ON` to cmake.
- The [solver/src/mpi](solver/src/mpi) directory a distributed memory or an hybrid
parallelisation version of the LBM algorithm. Activate it by passing the flags `-D_MPI=ON`
or `-D_MPI_OMP` to cmake.
- The [solver/src/cuda](solver/src/cuda) directory contains a GPU parallelisation version
of the LBM algorithm. Activate it by passing the flag `-D_CUDA=ON` to cmake.
- The [solver/src/mpi_cuda](solver/src/cuda) directory contains a multi-GPU
parallelisation version of the LBM algorithm. Activate it by passing the flag
`-D_CUDA_MPI=ON` to cmake.

Input parameters files are read from the [input](input) directory.
- The [input/inputJSON](input/inputJSON) directory contains the json files with
simulation parameters.
- The [input/inputPy](input/inputPy) directory contains the python scripts required
to convert the input files.

Output files are dumped in the [output](output) directory.
- The [output/outputVTR](output/outputVTR) directory contains all generated VTK
files.
- The [output/outputBackup](output/outputBackup) directory contains all generated
backup VTR files used to restart a simulation.

## Building metaLBM

Provided a working installation of CMake(>= 3.2) is present, `cd` to project
root and setup build directory

```shell
mkdir solver/build
cd solver/build
cmake ..
```

If the system's compiler is too old and does not support `C++11`, you can specify
a custom compiler:

```shell
cmake .. -DCMAKE_CXX_COMPILER=/path/to/compiler
```

To build the distributed-memory parallel version of metaLBM, a `MPI 3.0` compliant
library implementation is needed (eg `OpenMPI`, `MPICH`). If the library is already
installed, `CMake ` usually finds the library without problems.
In case multiple `MPI` libraries are installed, do

```shell
man mpicxx
```

to check if it is pointing to the right library. `CMake` uses it to determine
the correct library locations.

metaLBM uses extensively the [Boost](http://www.boost.org/) library and therefore a
working installation is required. You may have to set the environment variables
`BOOST_ROOT` to the location of the `Boost`'s root directory, `BOOST_INCLUDEDIR` to
the location of the `Boost`'s include directory, and `BOOST_LIBRARYDIR` to the
location of the `Boost`'s lib directory before running `cmake`:

```shell
export BOOST_ROOT="path/to/boost/root"
export BOOST_INCLUDEDIR="path/to/boost/include"
export BOOST_LIBRARYDIR="path/to/boost/lib"
```

In particular, metaLBM exploits the synergy between Boost Test and Boost Log to run
unit tests. To build all tests do:

```shell
make tests
```

If you want to add unit tests, just add a source file in `solver/test/`, edit the
correspondng `CMakeLists.txt` files to include your new targets and rebuild the
project `cmake ..` from your `solver/build` directory. There are predeclared functions
which make it easy to add new targets.

metaLBM supports scale-dependent forces, which are based on [FFTW](http://www.fftw.org)
library and therefore you will need a working installation. You may have to set the
environment variables `FFTW_INCLUDE_DIR` to the location of the `FFTW`'s include directory,
and `FFTW_LIBRARY_DIR` to the location of the `FFTW`'s lib directory before running `cmake`:

```shell
export FFTW_INCLUDE_DIR="path/to/fftw/include"
export FFTW_LIBRARY_DIR="path/to/fftw/lib"
```

To generate documentation for the project, you need a working `doxygen`
installation. Run from your `solverbuild` directory:

```shell
make doc
```

Documentation will be generated in a subfolder `solver/doc` and the webpage can be
found in `solver/doc/html/index.html`.

Check the [Wiki](https://gitlab.com/rooknrowl/metaLBM/wikis/home) for
additional tips on how to run a simulation.


## License
Please see [LICENSE](LICENSE).


## Contributing
Please see [CONTRIBUTING.md](CONTRIBUTING.md)
