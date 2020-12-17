CXX=nvcc
CXXFLAGS= -O3 -arch=compute_70

CXXFLAGS+= -std=c++14 -DPROTO_CUDA #-DPROTO_HIP #--expt-extended-lambda
CXXFLAGS+= -DsuperDebug

INCL+= -I../../include/ -I$(MVAPICH2_DIR)/include

SRCS := $(wildcard *.cu) 
OBJS=$(SRCS:.cu=)
all: $(OBJS)

$(OBJS): %: %.cu
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@.exe


clean:
	rm -rf *.exe *.o
