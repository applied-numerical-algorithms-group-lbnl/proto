CXX=nvcc
CXXFLAGS= -O3 -arch=compute_70

#CXXFLAGS+= -std=c++14 -DDIM=2 -lnvToolsExt #-DPROTO_HIP #--expt-extended-lambda
CXXFLAGS+= -g -std=c++11 -DPROTO_CUDA -DDIM=2 -x cu --expt-extended-lambda -lnvToolsExt 
#CXXFLAGS+= -DsuperDebug

INCL+= -I../../include/

SRCS := $(wildcard *.cu) 
OBJS=$(SRCS:.cu=)
all: $(OBJS)

$(OBJS): %: %.cu
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@.exe


clean:
	rm -rf *.exe *.o
