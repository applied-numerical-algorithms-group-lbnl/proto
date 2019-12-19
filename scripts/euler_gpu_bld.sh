#!/bin/bash
CC=nvcc
$CC -arch=compute_61 -O3 -g -DDIM=3 -DNUMCELLS=64 -DOMP_ENABLE=1 -DCUDA_ENABLE=1 -DDATAFLOW_ON=0 -DPROTO_CUDA=1 -std=c++11 -x cu -I./out -I./include -I./test -I./pdfl/src -I./examples/Euler/src -lomp examples/Euler/exec/Euler.cpp examples/Euler/src/EulerRK4.cpp examples/Euler/src/EulerOp.cpp -o euler_proto_gpu_n64.x
$CC -arch=compute_61 -O3 -g -DDIM=3 -DNUMCELLS=32 -DOMP_ENABLE=1 -DCUDA_ENABLE=1 -DDATAFLOW_ON=0 -DPROTO_CUDA=1 -std=c++11 -x cu -I./out -I./include -I./test -I./pdfl/src -I./examples/Euler/src -lomp examples/Euler/exec/Euler.cpp examples/Euler/src/EulerRK4.cpp examples/Euler/src/EulerOp.cpp -o euler_proto_gpu_n32.x
$CC -arch=compute_61 -O3 -g -DDIM=3 -DNUMCELLS=16 -DOMP_ENABLE=1 -DCUDA_ENABLE=1 -DDATAFLOW_ON=0 -DPROTO_CUDA=1 -std=c++11 -x cu -I./out -I./include -I./test -I./pdfl/src -I./examples/Euler/src -lomp examples/Euler/exec/Euler.cpp examples/Euler/src/EulerRK4.cpp examples/Euler/src/EulerOp.cpp -o euler_proto_gpu_n16.x

CC=pgc++
cp out/euler_step_3d_ser_acc.h out/euler_step.h
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=64 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ser_acc_n64.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=32 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ser_acc_n32.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=16 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ser_acc_n16.x

cp out/euler_step_3d_pf_acc.h out/euler_step.h
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=64 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_pf_acc_n64.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=32 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_pf_acc_n32.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=16 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_pf_acc_n16.x

cp out/euler_step_3d_ff_acc.h out/euler_step.h
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=64 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ff_acc_n64.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=32 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ff_acc_n32.x
$CC -fast -Minfo=accel -ta=tesla:cc60,fastmath,managed,beta -O3 -g -DDIM=3 -DNUMBOX=40 -DNUMCELLS=16 -DTILESIZE=1 -DOMP_ENABLE=1 -DPGI_ENABLE=1 -DDATAFLOW_ON=0 -DDATAFLOW_CODE=1 -I./out -I./include -I./test -I./pdfl/src -lomp ./test/euler_par.cpp -o euler_ff_acc_n16.x
