## Description

metaLBM is a C++11 header-only template Computational Fluid Dynamic software
based on the Lattice Boltzmann Method (LBM) that is addapted to both multi-CPUs and multi-GPUs architectures. 
metaLBM is meant to support multiple lattices and hile it can be used to perform various kind
of flow simulations, it has been developped with the aim of exploring turbulent flow
simulation and turbulence modelling within LBM.The library is hybrid parallel with
MPI for distributed memory and OpenMP threads for shared memory on multi-CPUs. It has 
been ported on multi-GPUs using CUDA and using the new NVSHMEM library for P2P communications 
(in collaboration with NVIDIA). It is under active development as of October 2018.

The [Lattice Boltzmann Method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods)
is a meso-scale approach to the simulation of fluid dynamics. Instead of solving the
Navier–Stokes equations, the discrete Boltzmann equation is solved to simulate the flow
of a Newtonian fluid with collision models such as Bhatnagar–Gross–Krook (BGK). By
simulating streaming and collision processes across a limited number of particles, the
intrinsic particle interactions evince a microcosm of viscous flow behavior applicable
across the greater mass.

metaLBM is to its authors knowledge, the only LBM code in full synergy with a FFT library. 
As turbulence is a multi-scale phenomenon, this allows to force, filter, and conduct analysis 
in Fourier space, which is based on [FFTW](http://www.fftw.org) library and therefore you will 
need a working installation.

The project provides a general framework for experimentation with simulations
of turbulent flows. It has an easy mechanism allowing to use various forcing schemes,
forces, and boundary conditions. The code currently runs on several medium size clusters
Stromboli at Bergische Universitaet Wuppertal (BUW) and Newturb at Università degli Studi
di Roma "Tor Vergata" and on larger supercomputers such as Galileo at CINECA in Bologna.
It is expected to run soon on JURECA and JURON at the Juelich Supercomputing Center (JSC).

## Overview

As of October 2018, metaLBM support the following features:

**Supported lattices**
- D1Q3
- D2Q5
- D2Q9
- D2Q13
- D2Q17
- D2Q21
- D3Q15
- D3Q19
- D3Q27
- D3Q33
- ... potentially any rectangular-shaped lattice

**Supported LBM algorithm**
- Fused collide and stream pull

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

**Supported forcing schemes**
- Shan-Chen
- Guo
- Exact-Difference Method

**Supported forces**
- Constant
- Sinusoidal
- Kolmogorov
- Spectral forcing on a shell of wavenumbers with time-independent phase (upcoming)
- Spectral forcing on a shell of wavenumbers with time-dependent phase (upcoming)

**Supported boundary conditions**
- Periodic BC
- Half-way bounceback (upcoming)
- Entropic bounceback (upcoming)

**Supported I/O**
- Paralell HDF5 

**Supported On-line analysis**
- Total energy, enstrophy, mach number
- Energy and forcing spectra
- Performances observables
- 
## Project Organisation

The project structure follows a common C++ layout within the [solver](solver) directory.

- The [test](test) directory contains source code for unit tests.
- The [include](include) directory contains headers-only library dependencies.
- The [src](src) directory contains all source files.
- The [doc](doc) directory contains documentation generating code
with doxygen.


## Building metaLBM

Building metaLBM requires a working installation of CMake(>= 3.9) (upcomming description).


To build tests, you will need a working intallation of
[GoogleTest](https://github.com/google/googletest) (upcomming).

If you want to add unit tests, just add a source file in [test/](test/),
edit the correspondng `CMakeLists.txt` files to include your new targets and rebuild the
project `cmake ..` from your [build](build) directory. There are
predeclared functions which make it easy to add new targets.

To generate documentation for the project, you need a working `doxygen`
installation.

Documentation will be generated in a subfolder [doc](doc) and the webpage
can be found in [doc/html/index.html](doc/html/index.html) (upcomming).

Check the [Wiki](https://gitlab.com/gtauzin/metaLBM/wikis/home) for
additional tips on how to build, compile, and run a simulation (upcomming).


## License
Please see [LICENSE](LICENSE).


## Contributing
Please see [CONTRIBUTING.md](CONTRIBUTING.md)
