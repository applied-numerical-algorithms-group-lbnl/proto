VARS_OLD := $(.VARIABLES)
#DEPXX=clang++
#CXX=clang++
DEPXX=mpicxx
CXX=mpicxx
DIM=3
USE_HDF5=TRUE
PROTO=../../..


LAPACK_LIB=/opt/homebrew/opt/openblas/lib
LAPACK_INC=/opt/homebrew/opt/openblas/include
#CXXOPTFLAGS=-g
CXXOPTFLAGS=-O2 -g
CXXOPTFLAGS+=-DmemoryFootPrint
#CXXOPTFLAGS+=-DPR_AMR
CXXOPTFLAGS+=-DPR_MMB
CXXOPTFLAGS+=-DPR_OPS
CXXOPTFLAGS+= -DPR_HDF5=TRUE
HDF5_LIB=/opt/homebrew/opt/hdf5-mpi/lib
HDF5_INC=/opt/homebrew/opt/hdf5-mpi/include
CXXOPTFLAGS+= -DPR_MPI=TRUE
CXXOPTFLAGS+= -DPR_STACK_ALLOC=4294967296
TARGET:=cubedSphereTest.exe

SRC=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp)) $(wildcard *.cpp) $(wildcard ../../_common/*.cpp) 
NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  ../../_common


CPPFLAGS:=-DDIM=$(DIM) -I$(PROTO)/include -I$(PROTO)/include/base -I$(PROTO)/include/amr -I$(PROTO)/include/mmb -I$(PROTO)/include/ops -I$(PROTO)/include/utils -I$(PROTO)/tests -I. -I../src -I../../_common
CXXFLAGS:= $(CXXOPTFLAGS) -std=c++20

ifdef HDF5_LIB
LIBS+= -L$(HDF5_LIB) -lhdf5 -lz
CPPFLAGS+= -I$(HDF5_INC)
endif

ifdef LAPACK_LIB
LIBS+= -L$(LAPACK_LIB) -lopenblas
CPPFLAGS+= -I$(LAPACK_INC)
endif

ifdef OLCF_HDF5_ROOT
LIBS+= -L$(OLCF_HDF5_ROOT)/lib -lhdf5 -lz
CPPFLAGS+= -I$(OLCF_HDF5_ROOT)/include
endif

ifdef FFTW3_DIR
LIBS+= -L$(FFTW3_DIR)/lib -lfftw3
CPPFLAGS+= -I$(FFTW3_DIR)/include
endif

DEP=$(DEPXX) -MM $(CXXFLAGS)
ifeq ($(DEPXX),nvcc)
DEP=$(DEPXX) -M $(CXXFLAGS)
LIBS+= -lnvToolsExt
endif

ifeq ($(DEPXX),hipcc)
CPPFLAGS+=-I$(ROCM_PATH)/rocthrust/include/ -I$(ROCM_PATH)/rocprim/include/ -I$(ROCM_PATH)/include/roctracer -I${MPICH_DIR}/include
CXXFLAGS+= -Wno-non-c-typedef-for-linkage -Wno-unused-command-line-argument
LIBS+= -L/$(ROCM_PATH)/lib/ -L${MPICH_DIR}/lib -lmpi -L${CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa
LIBS+= -lroctx64
endif


LINK=$(CXX)
OBJDIR = o
DEPDIR = d


$(TARGET): $(OBJ)
	$(LINK) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJDIR):
	@if [ ! -d "./$(OBJDIR)" ];then     \
		mkdir $(OBJDIR);           \
	fi

$(DEPDIR):
	@if [ ! -d "./$(DEPDIR)" ];then     \
		mkdir $(DEPDIR);           \
	fi

#compile non-main sources
$(OBJDIR)/%.o : %.cpp   $(OBJDIR) $(DEPDIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(CPPFLAGS)
	$(DEP) $<   $(CPPFLAGS) > $(DEPDIR)/$*.d



.PHONY: clean dataclean realclean print dirs


print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

clean:
	rm -rf d o *.o *.exe *.d ../../_common/*.o ../src/*.o ../../common/*.d ../src/*.d

dataclean:
	rm -rf *.hdf5 *.vtk *.table pout*



