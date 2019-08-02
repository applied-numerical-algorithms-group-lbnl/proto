#!/bin/bash

# Proto
cc  -DNUMCELLS=64 -DDIM=3                   -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp ./examples/Euler/src/EulerOp.cpp -o euler_step_3d.x -qopt-report-file=euler_step_3d.rpt

# Serial
cp out/euler_step_3d_ser_alloc.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_ser.x -qopt-report-file=euler_step_3d_ser.rpt
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_ser_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

# Serial-Tile
cp out/euler_step_3d_ser_tile_alloc.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_ser_tile.x -qopt-report-file=euler_step_3d_ser_tile.rpt
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_ser_tile_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

# Fuse
cp out/euler_step_3d_fuse_alloc.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_fuse.x -qopt-report-file=euler_step_3d_fuse.rpt
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_fuse_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

# Fuse-Tile
cp out/euler_step_3d_fuse_tile_alloc.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_fuse_tile.x -qopt-report-file=euler_step_3d_fuse_tile.rpt  
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_fuse_tile_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

# Fuse-Float
cp out/euler_step_3d_fuse_float.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_fuse_flt.x -qopt-report-file=euler_step_3d_fuse_flt.rpt
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_fuse_flt_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

# Dev (Fuse-All)
cp out/euler_step_3d_dev.h out/euler_step.h
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -g -dynamic -qopt-report=5 -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout ./test/euler_par.cpp -o euler_step_3d_dev.x -qopt-report-file=euler_step_3d_dev.rpt  
cc -DDIM=3 -DNUMCELLS=64 -DDATAFLOW_CODE=1 -DMPI_ENABLE=1 -DVTUNE_PERF=1 -g -dynamic -qopenmp -restrict -O3 -Iinclude -Iexamples/Euler/src -Ipdfl/src -Iout -o euler_step_3d_dev_vt.x ./test/euler_par.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify
