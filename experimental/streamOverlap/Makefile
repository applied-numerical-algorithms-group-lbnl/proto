CXX=nvcc
CXXFLAGS= -O3
CXXFLAGS+= -std=c++11 -arch=sm_70
INCL= -I${CUDA_HOME}/include -I{CUDA_HOME}/samples/common/inc
LIBS=

INC_MPI=-I/usr/common/software/mvapich2/2.3.2/pgi/19.10/cuda/10.2.89/include/
LIB_MPI=-L/usr/common/software/mvapich2/2.3.2/pgi/19.10/cuda/10.2.89/lib -lmpi


INCL+= $(INC_MPI)
LIBS+= $(LIB_MPI)

SRCS := $(wildcard *.cu) 
OBJS=$(SRCS:.cu=)
all: $(OBJS)

$(OBJS): %: %.cu
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@


clean:
	rm -rf main *.o
