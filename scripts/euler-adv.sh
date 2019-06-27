#!/bin/bash
if [ $SLURM_PROCID -eq 0 ]
then
 advixe-cl --collect $1 -exclude-files=Lists.hpp -flop -app-working-dir=$SCRATCH/proto --project-dir $SCRATCH/proto -- $2
else
  $1 
fi

