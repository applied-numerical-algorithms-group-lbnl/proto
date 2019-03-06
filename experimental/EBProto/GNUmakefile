
DEPXX=clang++
CXX=clang++
DIM=3
PROTO=../..

CXXOPTFLAGS=-g 
#CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+= --std=c++11
CXXOPTFLAGS+= -DDIM=$(DIM) -I$(PROTO)/include 

all : EBProto.exe

EBProto.exe : EBProto.cpp
	$(CXX) $(CXXOPTFLAGS) EBProto.cpp -o EBProto.exe

clean: 	
	rm *.exe
