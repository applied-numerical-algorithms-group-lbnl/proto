CXX=hipcc
CXXFLAGS= -Ofast

CXXFLAGS+= -Wnon-c-typedef-for-linkage -DDIM=3 -I$(PROTO)/include -I$(ROCM_PATH)/rocthrust/include/ -I$(ROCM_PATH)/rocprim/include/ 
CXXFLAGS+= -std=c++17 --amdgpu-target=gfx906 -DsuperDebug -DPROTO_CUDA -DPROTO_HIP #--expt-extended-lambda
LIBS=

INCL+= -I../../include/
LIBS+= $(LIB_MPI)

SRCS := $(wildcard *.cpp) 
OBJS=$(SRCS:.cpp=)
all: $(OBJS)

$(OBJS): %: %.cpp
	$(CXX) $(CXXFLAGS) $(INCL) $(LIBS) $< -o $@.exe


clean:
	rm -rf *.exe *.o
