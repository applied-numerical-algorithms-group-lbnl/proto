#!/bin/csh -f
if ($#argv != 4) then
  echo "Usage: kernel.sh dirname min_nx max_nx num_applies"
  exit 1
endif
set dirname = $1
set min_nx  = $2
set max_nx  = $3
set num_app = $4

echo "max_nx = $max_nx, min_nx = $min_nx, dirname = $dirname"

if (! -e $dirname) then
    set command = "mkdir $dirname"
    echo $command
    $command
endif

set nx_cur = $min_nx
while ($nx_cur <= $max_nx)
  set outfile = "$dirname/kernel$nx_cur.out"
  set command = "./applyKernel.exe -n $nx_cur -m $num_app >&  $outfile"
  echo $command
  $command
  set command = "mv proto.time.table $dirname/apply.nx.$nx_cur.time.table"
  echo $command
  $command

  set nx_cur = `expr 2 \* $nx_cur`
end

exit 0
