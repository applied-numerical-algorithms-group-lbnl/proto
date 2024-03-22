# Base Tools
Base contains all of the tools necessary for representing data on a union of rectangular patches on a single level as well as interfaces for defining pointwise functions, stencil operators, interpolation operators, and generalized finite volume operators.

## Classes
* Memory - An interface for platform independent memory manipulation.
* Point - An element of Z^DIM as an integer valued vector. Usable both on a CPU and GPU.
* Array - Built in version of std::array that can be used both on the CPU and GPU.
* Box - A logically rectangular region of Z^DIM. Does not maintain information on data centering.
* CoordPermutation - Represents a discrete rotation of coordinates
* Reduction - A tool for executing reduction operators on distributed datasets.
* BoxData - Representation of a single patch of data allocated on CPU or GPU. Very similar to FArrayBox in Chombo.
* DisjointBoxLayout - Representation of a collection of fixed-size Boxes which define a single level domain in Z^DIM. Contains tools for load balancing.
* LevelBoxData - A level worth of distributed data whose domain is defined by a DisjointBoxLayout. Can be allocted on CPU or GPU. Contains tools for copying and filling ghost cells.
* LevelOp - Extends the operability of the BoxOp interface to the scope of LevelData. Includes tools for applying boundary conditions through LevelBC (WIP)
* Stencil - Object representation of a stencil operation implemented as a sum of shift-coefficient pairs. Operates on the CPU or GPU.
* InterpStencil - An abstraction for an array of Stencils that are used to apply operators between datasets with non-trival refinement ratios (e.g. coarse-fine interpolation)
* Forall - A framework for defining and applying pointwise operator kernels on the CPU or the GPU.
* Operator - A collection of various namespace scoped utility functions for commonly used subroutines (e.g. convolution, product rule, etc.)
* HDF5Handler - Interface for I/O using HDF5. Contains tools for AMR and MMB libraries if they are used.

## Interfaces
* BoxOp - An extensible interface for defining generalized finite volume operators on a single patch of data.
* LevelBC - An extensible interface for defining boundary conditions (WIP).

