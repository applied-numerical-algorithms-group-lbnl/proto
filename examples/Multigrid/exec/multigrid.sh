#!/bin/csh
./Multigrid.exe -n 64 -r 1 -m 20 | tee  multigrid.out
mv proto.time.table _multigrid_bench/time.table.64
mv multigrid.out    _multigrid_bench/screen.out.64
./Multigrid.exe -n 128 -r 1 -m 20 | tee  multigrid.out
mv proto.time.table _multigrid_bench/time.table.128
mv multigrid.out    _multigrid_bench/screen.out.128
./Multigrid.exe -n 256 -r 1 -m 20 | tee  multigrid.out
mv proto.time.table _multigrid_bench/time.table.256
mv multigrid.out    _multigrid_bench/screen.out.256
./Multigrid.exe -n 512 -r 1 -m 20 | tee  multigrid.out
mv proto.time.table _multigrid_bench/time.table.512
mv multigrid.out    _multigrid_bench/screen.out.512

