
##DEPXX=clang++
##CXX=clang++
CXX=nvcc

DIM=2
PROTO=../..

#CXXOPTFLAGS=-g 
CXXOPTFLAGS=-g -G -DTHRUST_DEBUG -DPROTO_CUDA -x cu
CXXOPTFLAGS+= --std=c++11
CXXOPTFLAGS+= -DDIM=$(DIM) -I$(PROTO)/include -I$(PROTO)/EBProto/include 

all : func_mess.exe ebforall.exe

func_mess.exe : func_mess.cpp 
	$(CXX) $(CXXOPTFLAGS) func_mess.cpp -o func_mess.exe

ebforall.exe : EBForAll.cpp 
	$(CXX) $(CXXOPTFLAGS) EBForAll.cpp -o ebforall.exe

clean: 	
	rm -rf *.exe *.dSYM
