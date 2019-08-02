#!/bin/bash
#SBATCH --partition=debug
#SBATCH --perf=vtune
#SBATCH --constraint="haswell"
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=euler_run
#SBATCH -o euler_run.%j       # output and error file name (%j expands to jobID)
#SBATCH --mail-user=edavis@lbl.gov
#SBATCH --mail-type=ALL

n=64
t=1

SDE='sde64 -hsw'

module load sde/8.35
module load vtune/2018.up2

export KMP_AFFINITY=granularity=core,compact,1
export OMP_NUM_THREADS=$t

# Normal runs
srun -n 64 ./euler_step_3d.x 7
srun -n 64 ./euler_step_3d_ser.x 7
srun -n 64 ./euler_step_3d_ser_tile.x 7
srun -n 64 ./euler_step_3d_fuse.x 7
srun -n 64 ./euler_step_3d_fuse_tile.x 7
srun -n 64 ./euler_step_3d_fuse_flt.x 7
srun -n 64 ./euler_step_3d_dev.x 7

# SDE Runs
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_ser -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_ser_vt.x 1
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_ser_tile -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_ser_tile_vt.x 1
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_fuse -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_fuse_vt.x 1
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_fuse_tile -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_fuse_tile_vt.x 1
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_fuse_flt -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_fuse_flt_vt.x 1
srun -n 64 -c 1 sde64 -hsw -d -iform 1 -omix euler_step_3d_dev -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- ./euler_step_3d_dev_vt.x 1

# VTune runs
srun -n 64 -c 1 amplxe-cl -start-paused -r euler_step_3d_ser -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi -- ./euler_step_3d_ser_vt.x 1
srun -n 64 -c 1 amplxe-cl -start-paused -r euler_step_3d_ser_tile -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi -- ./euler_step_3d_ser_tile_vt.x 1
srun -n 64 -c 1 amplxe-cl -start-paused -r euler_step_3d_fuse -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi -- ./euler_step_3d_fuse_vt.x 1
srun -n 64 -c 1 amplxe-cl -start-paused -r euler_step_3d_fuse_tile -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi -- ./euler_step_3d_fuse_tile_vt.x 1
srun -n 64 -c 1 amplxe-cl -start-paused -r euler_step_3d_dev -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi -- ./euler_step_3d_dev_vt.x 1
