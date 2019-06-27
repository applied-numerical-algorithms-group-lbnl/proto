#!/bin/bash
cc -DVTUNEPERF -Ipdfl/src -Iout -g -dynamic -qopenmp -O3 -o euler test/eulerpar.cpp -I$VTUNE_AMPLIFIER_XE_2018_DIR/include -L$VTUNE_AMPLIFIER_XE_2018_DIR/lib64 -littnotify

