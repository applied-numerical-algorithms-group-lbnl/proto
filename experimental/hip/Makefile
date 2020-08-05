CXX=hipcc
CXXFLAGS= -Ofast

CXXFLAGS+= -Wnon-c-typedef-for-linkage -Wno-unused-result
CXXFLAGS+= -std=c++17 --amdgpu-target=gfx906 -DPROTO_CUDA -DPROTO_HIP #--expt-extended-lambda
#CXXFLAGS+= -DsuperDebug

INCL+= -I../../include/ -I$(MPI_HOME)/include
LIBS+= -L$(MPI_HOME)/lib -lmpi

SRCS := $(wildcard *.cu) 
OBJS=$(SRCS:.cu=)
all: $(OBJS)

$(OBJS): %: %.cu
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@.exe


clean:
	rm -rf *.exe *.o
