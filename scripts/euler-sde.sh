#!/bin/bash
#SBATCH --partition=debug
#SBATCH --perf=vtune
#SBATCH --constraint="haswell"
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --job-name=euler_fuse
#SBATCH -o euler_fuse.%j       # output and error file name (%j expands to jobID)
#SBATCH --mail-user=edavis@lbl.gov
#SBATCH --mail-type=ALL

n=64
t=1

SDE='sde64 -hsw'

module load sde/8.35
module load vtune/2018.up2

echo $SDE
export CH_TIMER=1
baseline=yes
sde=yes
vtune=yes
export KMP_AFFINITY=granularity=core,compact,1
exe=./euler_fuse
suffix=${t}t_${n}n
export OMP_NUM_THREADS=$t

if [ "$sde" == "yes" ]; then
  echo ""
  echo "--------------------------------------------------"
  echo "----->> Running w/SDE <<-----"
  echo "--------------------------------------------------"
  srun -n $n -c $t $SDE -d -iform 1 -omix sde_${suffix}.out -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- $exe
  echo "----->> Generating SDE Report <<-----"
  echo "For performance, the SDE report is best done on an external login node"
  echo "Run the following command: "
  echo "\$ ./parse-sde.sh sde_${suffix}.out*"
fi

if [ "$vtune" == "yes" ]; then
  echo ""
  echo "--------------------------------------------------"
  echo "----->> Running w/Vtune <<-----"
  echo "--------------------------------------------------"
  srun -n $n -c $t amplxe-cl -start-paused -r vtbw_${suffix} -data-limit=0 -collect memory-access -finalization-mode=none -trace-mpi --  $exe
  echo "----->> Finalizing VTune and generating report <<-----"
  echo "For performance, the finalize and report are best done on an external login node"
  echo "Run the following commands: "
  echo "\$ amplxe-cl -report summary -r vtbw_${suffix} > vtbw_${suffix}.summary"
  echo "\$ ./parse-vtune.shl vtbw_${suffix}.summary"
  echo "For a report breaking out individual DRAM channels per socket (packages):"
  echo "\$ amplxe-cl -report hw-events -r vtbw_${suffix}"
fi

