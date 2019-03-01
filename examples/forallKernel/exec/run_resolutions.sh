#!/bin/bash -f
if [ "$#" -ne 4 ] 
then
  echo "Usage: kernel.sh dirname min_nx max_nx num_applies"
  exit 1
fi
dirname=$1
min_nx=$2
max_nx=$3
num_app=$4

echo "max_nx = $max_nx, min_nx = $min_nx, dirname = $dirname"

if [ ! -d $dirname ] 
then
    command="mkdir $dirname"
    echo $command
    $command
fi

nx_cur=$min_nx
while [ $nx_cur -le $max_nx ]
do
  outfile="$dirname/kernel$nx_cur.out"
  command="./forallKernel.exe -n $nx_cur -m $num_app >&  $outfile"
  echo $command
  $command
  command="mv proto.time.table $dirname/forall.nx.$nx_cur.time.table"
  echo $command
  $command

  nx_cur=`expr 2 \* $nx_cur`
done
