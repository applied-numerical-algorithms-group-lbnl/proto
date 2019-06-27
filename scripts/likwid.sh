#!/bin/bash
salloc -N 1 -q interactive -C haswell -t 15:00 --perf=likwid
module load likwid
cc -DLIKWID_PERFMON -I$LIKWID_INCLUDE -I./pdfl/src -Iout -L$LIKWID_LIB -llikwid -dynamic ./test/eulerpar.cpp -o ./euler
likwid-mpirun -pin 0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15 -g FLOPS_AVX -m -O ./euler_ser > euler_ser-likwid-avx.csv

