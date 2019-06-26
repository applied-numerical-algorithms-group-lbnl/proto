VARS_OLD := $(.VARIABLES)
DEPXX=nvcc
CXX=nvcc
#DEPXX=g++
#CXX=g++
DIM=3
PROTO=../..

#CXXOPTFLAGS=-g
CXXOPTFLAGS=-O3 
CXXOPTFLAGS+=-x cu --expt-extended-lambda
CXXOPTFLAGS+=-DPROTO_CUDA=1
TARGET:=streams.exe


SRC=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp))  $(wildcard *.cpp) $(wildcard ../../common/*.cpp) 
OBJ=$(subst .cpp,.o, $(SRC))

DEPENDS:=$(subst .cpp,.d,$(SRC))
VPATH=. ../src  ../../../common


CPPFLAGS:=-DDIM=$(DIM) -I$(PROTO)/include  -I. -I../src -I../../common
CXXFLAGS:= $(CXXOPTFLAGS) -std=c++11

ifdef FFTW3_DIR
LIBS+= -L$(FFTW3_DIR)/lib -lfftw3
CPPFLAGS+= -I$(FFTW3_DIR)/include
endif

DEP=$(DEPXX) -MM $(CXXFLAGS)
ifeq ($(DEPXX),nvcc)
DEP=$(DEPXX) -M $(CXXFLAGS)
LIBS+= -lnvToolsExt
endif

LINK=$(CXX)

$(TARGET) : $(OBJ)
	$(LINK) -o $(TARGET) $(OBJ) $(LIBS)

#compile non-main sources
%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS)
	$(DEP) $<   $(CPPFLAGS) > $*.d

.PHONY: clean realclean print dirs

clean:
	rm -f *.o *.exe *.d ../../common/*.o $(EULERSRC)/*.o

print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

#
#DEPXX=nvcc
#CXX=nvcc
##DEPXX=clang++
##CXX=clang++
#DIM=3
#PROTO=../..
#
#CXXOPTFLAGS=-O3 
##CXXOPTFLAGS=-g -G -DTHRUST_DEBUG
#CXXOPTFLAGS+= --std=c++11 
#CXXOPTFLAGS+=-DPROTO_CUDA=1 -DDIM=$(DIM) -x cu --expt-extended-lambda -I$(PROTO)/include 
#
#all : streams.exe
#
#streams.exe : streams.o EulerOp.o
#	$(CXX) streams.o EulerOp.o -o streams.exe
#
#streams.o : streams.cpp
#	$(CXX) $(CXXOPTFLAGS) -c stream.cpp -o stream.o
#
#EulerOp.o : ../../examples/Euler/src/EulerOp.cpp ../../examples/Euler/src/EulerOp.H
#	$(CXX) $(CXXOPTFLAGS) -I../../examples/Euler/src -c ../../examples/Euler/src/EulerOp.cpp
#
#clean: 	
#	rm *.exe
