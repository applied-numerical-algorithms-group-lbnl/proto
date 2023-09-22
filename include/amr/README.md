# AMR Tools
Contains all of the tools necessary for representing data and operators on an AMR hierarchy.

## Classes
* AMRGrid - Abstraction for an array of DisjointBoxLayouts that define the domain of an AMR Hierarchy.
* LevelFluxRegister - Special data holder which holds data on coarse-fine boundaries of an AMRHierarchy and facilitates the refluxing operation
* AMRData - Abstraction for an array of LevelBoxData defined on an AMRGrid. Includes tools for regridding. 
* AMROp - Extends the operability of a BoxOp to the scope of an AMRData. Includes tools for coarse-fine boundary refluxing.
