#This is a standalone Makefile for this example because Proto configure doesn't handle
#compiling parallel applications yet.
#In general, the user only needs to set CXX for a given machine. The MPI compiler will
#compile and link with the proper MPI files. However, here are some per-machine instructions
#to make sure everything is set up correctly.
#1) ANAG machine (tried on Madonna)
#   a) Set CXX=mpicxx (Set DEPXX to CXX)
#   b) The ANAG machines usually load openmpi by default. Check that mpirun and mpicxx point to 
#      the same openmpi directory. Can change to mpich using "module unload openmpi; module load mpich"
#2) Cori (Intel programming environment + Haswell node)
#   a) The Intel programming environment is loaded by default. Check using "module list". Can switch
#      to GNU or Cray, but we haven't checked this makefile works for these environments. The proper
#      MPI modules should be loaded when you load a particular environment. Nothing else needs to be done.
#   b) Set CXX=CC (Set DEPXX to CXX). CC is Cori's wrapper for the underlying compilers. 
#3) Mac (ask Phil)
#   a) Set CXX=mpic++ and DEPXX=clang++ (this is the version Phil has compiled)

VARS_OLD := $(.VARIABLES)
PR_MPI = TRUE
ifdef PR_MPI
CXX=mpicxx
else
CXX=g++
endif
DEPXX=$(CXX)
DIM=3
PROTO=../../..

#CXXOPTFLAGS=-g 
CXXOPTFLAGS=-O3
TARGET:=EulerLevelData.exe
SRC:=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp)) $(wildcard ../../common/*.cpp)
SRC += EulerTest.cpp
NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  ../../common


CPPFLAGS:=-DDIM=$(DIM) -DPR_MPI=$(PR_MPI) -I$(PROTO)/include  -I$(PROTO)/EBProto/include -I. -I../src -I../../common
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
OBJDIR = o
DEPDIR = d


$(TARGET): $(OBJ)
	$(LINK) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJDIR):
	if [ ! -d "./$(OBJDIR)" ];then     \
		mkdir $(OBJDIR);           \
	fi

$(DEPDIR):
	if [ ! -d "./$(DEPDIR)" ];then     \
		mkdir $(DEPDIR);           \
	fi

#compile non-main sources
$(OBJDIR)/%.o : %.cpp   $(OBJDIR) $(DEPDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS)
	$(DEP) $<   $(CPPFLAGS) > $(DEPDIR)/$*.d



.PHONY: clean realclean print dirs


print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

clean:
	rm -rf d o *.o *.exe *.d ../../common/*.o ../src/*.o ../../common/*.d ../src/*.d




