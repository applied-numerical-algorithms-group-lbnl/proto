CXX=hipcc
CXXFLAGS= -Ofast

CXXFLAGS+= -Wnon-c-typedef-for-linkage
CXXFLAGS+= -std=c++17 --amdgpu-target=gfx906 -DsuperDebug -DPROTO_CUDA -DPROTO_HIP #--expt-extended-lambda
#CXXFLAGS+= -std=c++11 -v --amdgpu-target=gfx906 -DPROTO_CUDA -DPROTO_HIP #--expt-extended-lambda
LIBS=

INCL+= -I../../include/
LIBS+= $(LIB_MPI)

SRCS := $(wildcard *.cu) 
OBJS=$(SRCS:.cu=)
all: $(OBJS)

$(OBJS): %: %.cu
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@.exe


clean:
	rm -rf *.exe *.o
