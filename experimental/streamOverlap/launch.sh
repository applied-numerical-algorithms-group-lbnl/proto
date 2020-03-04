export OMP_NUM_THREADS=4

rm aa
make clean
make
srun nvprof -o aa ./main 
nvvp aa
