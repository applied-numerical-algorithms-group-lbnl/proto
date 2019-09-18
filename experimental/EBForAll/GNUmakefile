VARS_OLD := $(.VARIABLES)
DEPXX=nvcc
CXX=nvcc
DIM=2
PROTO=../..

CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
CXXOPTFLAGS+=-x cu --expt-extended-lambda -std=c++11
CXXOPTFLAGS+=-DPROTO_CUDA=1
TARGET:=iota.exe

CPPFLAGS:=-DDIM=$(DIM) -I$(PROTO)/include  -I$(PROTO)/EBProto/include


all : func_mess.exe ebforall.exe iota.exe

iota.exe : iota.cpp 
	$(CXX) $(CXXOPTFLAGS) $(CPPFLAGS) iota.cpp -o iota.exe

func_mess.exe : func_mess.cpp 
	$(CXX) $(CXXOPTFLAGS) $(CPPFLAGS) func_mess.cpp -o func_mess.exe

ebforall.exe : EBForAll.cpp 
	$(CXX) $(CXXOPTFLAGS) $(CPPFLAGS) EBForAll.cpp -o ebforall.exe

clean: 	
	rm -rf *.exe *.dSYM

