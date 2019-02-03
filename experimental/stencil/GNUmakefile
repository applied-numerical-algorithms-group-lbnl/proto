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

all : applyStream.exe forallStream.exe forallStruct.exe scaleCompare.exe forallZincStruct.exe doNothingForAll.exe forallOneD.exe baroqueStencil.exe forallHardWired.exe baroqueForAllMin.exe forallSimple.exe

applyStream.exe : applyStream.cu
	$(CXX) $(CXXOPTFLAGS) applyStream.cu -o applyStream.exe

forallStream.exe : forallStream.cu
	$(CXX) $(CXXOPTFLAGS) forallStream.cu -o forallStream.exe

forallStruct.exe : forallStruct.cu
	$(CXX) $(CXXOPTFLAGS) forallStruct.cu -o forallStruct.exe

scaleCompare.exe : scaleCompare.cu
	$(CXX) $(CXXOPTFLAGS) scaleCompare.cu -o scaleCompare.exe

forallZincStruct.exe : forallZincStruct.cu
	$(CXX) $(CXXOPTFLAGS) forallZincStruct.cu -o forallZincStruct.exe

doNothingForAll.exe : doNothingForAll.cu
	$(CXX) $(CXXOPTFLAGS) doNothingForAll.cu -o doNothingForAll.exe

forallOneD.exe : forallOneD.cu
	$(CXX) $(CXXOPTFLAGS) forallOneD.cu -o forallOneD.exe

baroqueStencil.exe : baroqueStencil.cu
	$(CXX) $(CXXOPTFLAGS) baroqueStencil.cu -o baroqueStencil.exe

forallHardWired.exe : forallHardWired.cu
	$(CXX) $(CXXOPTFLAGS) forallHardWired.cu -o forallHardWired.exe

baroqueForAllMin.exe : baroqueForAllMin.cu
	$(CXX) $(CXXOPTFLAGS) baroqueForAllMin.cu -o baroqueForAllMin.exe

forallSimple.exe : forallSimple.cu
	$(CXX) $(CXXOPTFLAGS) forallSimple.cu -o forallSimple.exe

clean: 	
	rm *.exe
