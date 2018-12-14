DEPXX=nvcc
CXX=nvcc
DIM=3
PROTO=../..

#CXXOPTFLAGS=-O3 
CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+=-x cu --expt-extended-lambda --std=c++11
CXXOPTFLAGS+=-DPROTO_CUDA=1 -DDIM=$(DIM) -I$(PROTO)/include 

all : applyStream.exe

applyStream.exe : applyStream.cu
	$(CXX) $(CXXOPTFLAGS) applyStream.cu -o applyStream.exe

clean: 	
	rm *.exe
