DEPXX=nvcc
CXX=nvcc
#DEPXX=clang++
#CXX=clang++
DIM=3
PROTO=../..

CXXOPTFLAGS=-O3 
#CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+= --std=c++11
CXXOPTFLAGS+=-DPROTO_CUDA=1 -DDIM=$(DIM) -I$(PROTO)/include 

all : applyStream.exe forallStream.exe

applyStream.exe : applyStream.cu
	$(CXX) $(CXXOPTFLAGS) applyStream.cu -o applyStream.exe

forallStream.exe : forallStream.cu
	$(CXX) $(CXXOPTFLAGS) forallStream.cu -o forallStream.exe

clean: 	
	rm *.exe
