# Proto
Courtesy of the Applied Numerical Algorithms Group, Lawrence Berkeley National Laboratory

### Introduction
Proto is a flexible, expressive, and lightweight performance library focused on providing tools for building high-order simulations on structured grids using adaptive mesh refinement (AMR). Proto provides:
* A clean, multi-level interface allowing developers to use as much or little of Proto as is needed for their application
* High level abstractions for easily building AMR-based solvers on simple or mapped-multiblock problem domains
* Expressive methods for specifying complex, high order computations using a minimalist, math-like syntax while avoiding excessive boilerplate
* A flexible backend equipped to target various heterogeneous computing platforms containing CPUs as well as both NVidia and AMD GPUs

## Including Proto
Proto requires very little effort to include in an external project. Proto consists entirely of header files, so there is no "build" at all. The only requirements are to #include "Proto.H" in the top level application, link the proto/include directory in the external application's build, and define any relevant compile time variables (see below). Do not #include any other headers from Proto other than Proto.H as doing so may lead to linking issues.
 ### Compile Time Constants
 * DIM - Maximum number of spatial dimensions
 * PR_MPI - If defined, MPI code will be toggled on. Compatible with either of the following GPU backend options
 * PROTO_CUDA - If defined, the cuda GPU compatible backend will be toggled on. Not compatible with PROTO_HIP
 * PROTO_HIP - If defined, the HIP GPU compatible backend will be toggled on. Not compatible with PROTO_CUDA
 * PR_HDF5 - If defined, HDF5 I/O tools will be included. 
 * PR_AMR - If defined, adaptive mesh refinement tools will be included
 * PR_MMB - If defined, mapped multiblock tools will be included. Requires PR_OPS
 * PR_OPS - If defined, built in linear algebra tools will be included. Requires LAPACK.
 * PR_BLIS - If defined, uses the third party BLIS library to handle linear algebra instead of LAPACK. Requires PR_OPS. (Interface is experimental)

## Build instructions:
The following instructions pertain to building applications Proto's native examples and test files. External applications using Proto as a third party library should follow the instructions in the previous section. 

### CMake version
* The minumum required CMake version is 3.20. If an older version is the default on your platform, load a `cmake` module with a sufficiently new version. Configuration has been tested through CMake 3.21. Users who are not comfortable with using CMake are invited to use the provided proto_make wrapper script to call CMake using a GNUMake style interface (see the proto_make section below).
### Cori Modules
* If doing a CUDA build, load the `cudatoolkit` module
* If doing an MPI build, load the `openmpi` module
* If using HDF5, load the `cray-hdf5` module
### Summit Modules
* If doing a CUDA build, load the `cuda/11` module. If you get configuration fails with a message like `error: identifier "__ieee128" is
  undefined` you may have an older version of CUDA loaded.
* If doing an MPI build, load the `spectrum-mpi` module
* If using HDF5, load the `hdf5` module

### Submodules
* This repository uses the [BLT](https://github.com/LLNL/blt) library of CMake macros during configuration, and contains it as a submodule.
* If you *haven't* yet cloned this repository, include the `--recurse-submodules` flag in your `git clone` command so that the submodule will be initialized and updated.
* If you *have* already cloned this repository, run these two commands to checkout the BLT commit that this project links to:
   - `git submodule init`
   - `git submodule update`
   
   After doing this, the `.gitmodules` file will show the path and URL of the submodule.

### Configuring
* If using CMake, the simplest command assumes you are at the top directory of the repo and is `cmake -B build`. More generically, the command is `cmake -S <source-dir> -B <build-dir>`. The argument to `-S` is the source directory containing the top-level `CMakeLists.txt` file: if not provided, it is assumed to be the current directory. The argument to `-B` is the directory where the binaries should be built, and does not need to exist when the `cmake` command is invoked. Additionally, there are various options which can be set during this step, each preceded by the `-D` flag. Valid choices are listed in brackets, with defaults in italics. They are:
   - ENABLE_CUDA=[ON, *OFF*]
   - ENABLE_HIP=[ON, *OFF*]
   - ENABLE_MPI=[ON, *OFF*]
   - ENABLE_HDF5=[*ON*, OFF]
   - ENABLE_OPENMP=[ON, *OFF*]
   - Build executables from the `examples` subdirectory: ENABLE_EXAMPLES=[*ON*, OFF]
   - Build default executables from the `tests` subdirectory: ENABLE_TESTS=[*ON*, OFF]
   - Build all from the `tests` subdirectory: ENABLE_ALL_TESTS=[ON, *OFF*]
   - Floating point precision: PREC=[SINGLE, *DOUBLE*]
   - Dimensionality of examples: DIM=[*2*, 3]
   - Optimization level: CMAKE_BUILD_TYPE=[*Debug*, Release, MinSizeRel, RelWithDebInfo]
   - Size of allocations from the stack: STACK=\<int\> (default is 4GB)
   - Turn on timers: TIMERS=[*ON*, OFF]
   - Turn on code that checks if copying/aliasing is working correctly: MEMCHECK=[ON, *OFF*]
   - Print the amount of data allocated per protoMalloc: MEMTRACK=[ON, *OFF*]

#### `proto_make`
* The `proto_make` Python script offers ease-of-use configuring for those unfamiliar with CMake. 
Using Python to invoke the script will display all targets within the Proto project -- 
some of these are examples while others are tests written for CI purposes. A target
for the script to compile can be specified using the `-t` flag. Special targets are `ALL` and `TEST`.
The former will build all targets within the project while the latter will take the additional step
of running all CI tests. If specific targets are provided (even if they are CI tests) they are only
compiled, not run. However, they have a symlink created for their executable in the root source directory
with an `.exe` suffix.

The full set of options available during CMake configuration are available when using `proto_make` too, 
but without any `ENABLE_` prefixes.
   
### Building
* To build all of the CMake targets in this project in parallel, run `cmake --build <build-dir> --parallel`. An integer can be provided at the end of this command to set the level of parallelization. The `<build-dir>` is NOT the name of the source directory containing the targets you want built (such as `examples`), but rather the name of the directory provided to the `-B` flag in the configuration command.
* `VERBOSE=1` may be added to the end of the build commend to view compilation details
* A specific target may be specified by adding `--target <target-name>` to the build command.
* The build command is simplified if you first navigate inside to build-dir. In this environment, the build command is simply `make <target-name>`
* After moving to the build directory, run all compiled tests by giving the command `make test` to a SLURM script or `salloc` line. This command must be invoked on the proper node (eg. one with GPUs if doing a CUDA run) for the tests to run successfully. The BLT submodule used by this project has various "smoke" tests that can be used to determine if the tests are being run on the proper architecture for the desired configuration. They are automatically included in the `test` target.
### Specific Architectures
| Macbooks with AMR processors | 
| ---------------------------- |
| An additional flag needs to be added to the configure command to specify the correct target architecture: |
| `-D CMAKE_APPLE_SILICON_PROCESSOR=arm64` |
| NERSC Cori GPU | 
| -------------- |
Modules: 
    - cgpu
    - cuda/11
Details:
From the regular Cori node load the `cgpu` module. Then, from a Cori GPU node (accessed either through a SLURM script or interactively) load the `cuda/11` module.
| NERSC Perlmutter | 
| ---------------- |
Modules:
    - cmake
    - cudatoolkit
    - cray-hdf5-parallel
    - craype-accel-nvidia80
Environment Variables:
    - `MPICH_GPU_SUPPORT_ENABLED=1`
| OLCF Crusher | 
| ------------ |
The following modules are needed when compiling on OLCF's Crusher machine:
    - cmake
    - hdf5
    - rocm
    - PrgEnv-amd
    - craype-accel-amd-gfx90a

### Cleanup
The script `cleanup.sh` can be run to remove output files in the root source generated when running the project's targets.
To remove output files from the build directory run `cmake -B <build-dir> -t dataclean`.
To remove executables from the build directory run `cmake -B <build-dir> -t clean`.
