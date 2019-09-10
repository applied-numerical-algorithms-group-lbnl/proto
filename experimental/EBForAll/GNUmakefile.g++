
##DEPXX=clang++
##CXX=clang++
##CXX=nvcc
CXX=g++
DIM=2
PROTO=../..

CXXOPTFLAGS=-g 
#CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+= --std=c++11
CXXOPTFLAGS+= -DDIM=$(DIM) -I$(PROTO)/include -I$(PROTO)/EBProto/include 

all : ebforall.exe

ebforall.exe : EBForAll.cpp 
	$(CXX) $(CXXOPTFLAGS) EBForAll.cpp -o ebforall.exe

clean: 	
	rm -rf *.exe *.dSYM
