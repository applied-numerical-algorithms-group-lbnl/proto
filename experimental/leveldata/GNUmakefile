DEPXX=nvcc
CXX=nvcc
#DEPXX=clang++
#CXX=clang++
DIM=3
PROTO=../..

CXXOPTFLAGS=-O3 
#CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+= --std=c++11 -ccbin clang
CXXOPTFLAGS+=-DPROTO_CUDA=1 -DDIM=$(DIM) -x cu --expt-extended-lambda -I$(PROTO)/include 

all : streams.exe

streams.exe : streams.o EulerOp.o
	$(CXX) streams.o EulerOp.o -o streams.exe

streams.o : streams.cpp
	$(CXX) $(CXXOPTFLAGS) -c stream.cpp -o stream.o

EulerOp.o : ../../examples/Euler/src/EulerOp.cpp ../../examples/Euler/src/EulerOp.H
	$(CXX) $(CXXOPTFLAGS) -I../../examples/Euler/src -c ../../examples/Euler/src/EulerOp.cpp

clean: 	
	rm *.exe
