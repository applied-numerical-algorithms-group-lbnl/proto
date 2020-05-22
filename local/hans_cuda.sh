##!/bin/bash

echo "Running script for Cori GPU build of proto/brick ..."

echo "Loading modules for CUDA build..."
module purge
module load esslurm
module load cmake
module load intel/18.0.1.163 cuda/10.1.168 mvapich2/2.3.1 hdf5-parallel/1.10.1-pgi

echo "./configure --dim 2 --opt OPT --timers TRUE --cxx nvcc"
./configure --dim 2 --opt OPT --timers TRUE --cxx nvcc

echo "getting bricklib:\n git submodule update --init"
git submodule update --init

echo "cmake in build_cuda ..."
mkdir build_cuda; cd build_cuda
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=ON
# cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=ON

echo "make laplacian"
make laplacian

echo "if you're missing mpi headers, add mvapich2 to brick flags.makes"

echo "export PR_TIME=TRUE"
export PR_TIME=TRUE
