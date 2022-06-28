# Proto
Courtesy of the Applied Numerical Algorithms Group, Lawrence Berkeley National Laboratory

### Introduction
Proto is a middleware layer that allows performance portability in scientific computations on structured grids. 
* Proto has a high-level interface to enable stencil computations and pointwise calculations with a minimum of code.  This productivity layer makes for much cleaner and smaller scientific codes.
* Proto allows users to have performance portability between platformas (including accelerator-based platforms) without having to change user code.
* Proto includes data structures which support embedded boundary (EB) calculations.  Embedded boundary calculations done properly require data to live on a more complex graph structure than simple arrays can support. 

## Build instructions:
### CMake version
* The minumum required CMake version is 3.20. If an older version is the default on your platform, load a `cmake` module with a sufficiently new version. Configuration has been tested through CMake 3.21.
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
* The simplest command assumes you are at the top directory of the repo and is `cmake -B build`. More generically, the command is `cmake -S <source-dir> -B <build-dir>`. The argument to `-S` is the source directory containing the top-level `CMakeLists.txt` file: if not provided, it is assumed to be the current directory. The argument to `-B` is the directory where the binaries should be built, and does not need to exist when the `cmake` command is invoked. Additionally, there are various options which can be set during this step, each preceded by the `-D` flag. Valid choices are listed in brackets, with defaults in italics. They are:
   - ENABLE_CUDA=[ON, *OFF*]
   - ENABLE_MPI=[ON, *OFF*]
   - ENABLE_HDF5=[*ON*, OFF]
   - ENABLE_OPENMP=[ON, *OFF*]
   - Build executables from the `examples` subdirectory: ENABLE_EXAMPLES=[*ON*, OFF]
   - Build executables from the `tests` subdirectory: ENABLE_TESTS=[*ON*, OFF]
   - Floating point precision: PREC=[SINGLE, *DOUBLE*]
   - Dimensionality of examples: DIM=[*2*, 3]
   - Optimization level: CMAKE_BUILD_TYPE=[*Debug*, Release, MinSizeRel, RelWithDebInfo]
   - Size of allocations from the stack: STACK=\<int\> (default is 4GB)
   - Turn on timers: TIMERS=[*ON*, OFF]
   - Turn on code that checks if copying/aliasing is working correctly: MEMCHECK=[ON, *OFF*]
   - Print the amount of data allocated per protoMalloc: MEMTRACK=[ON, *OFF*]

| For Macbooks with ARM processors |
| -------------------------------- |
| An additional flag needs to be added to the configure command to specify the correct target architecture: |
| `-DCMAKE_APPLE_SILICON_PROCESSOR=arm64` |
   
### Building
* To build all of the CMake targets in this project in parallel, run `cmake --build <build-dir> --parallel`. An integer can be provided at the end of this command to set the level of parallelization. The `<build-dir>` is NOT the name of the source directory containing the targets you want built (such as `examples`), but rather the name of the directory provided to the `-B` flag in the configuration command.
   
### Testing
| For Cori GPU only| 
| ---------------- |
| Extra modules are needed to run CUDA-enabled executables on Cori. From the regular Cori node load the `cgpu` module. Then, from a Cori GPU node (accessed either through a SLURM script or interactively) load the `cuda/11` module. |
* After moving to the build directory, run all compiled tests by giving the command `make test` to a SLURM script or `salloc` line. This command must be invoked on the proper node (eg. one with GPUs if doing a CUDA run) for the tests to run successfully. The BLT submodule used by this project has various "smoke" tests that can be used to determine if the tests are being run on the proper architecture for the desired configuration. They are automatically included in the `test` target.
