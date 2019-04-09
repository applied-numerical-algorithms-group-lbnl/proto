VARS_OLD := $(.VARIABLES)
DEPXX=g++
CXX=g++
DIM=2
PROTO=../../../..
EBPROTO=../../../../EBProto

CXXOPTFLAGS=-g
TARGET:=GeometryConvergence.exe


SRC=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp)) $(wildcard *.cpp) $(wildcard ../../common/*.cpp) 
OBJ=$(subst .cpp,.o, $(SRC))

DEPENDS:=$(subst .cpp,.d,$(SRC))
VPATH=. ../src  ../../../common


CPPFLAGS:=-DDIM=$(DIM) -I$(EBPROTO)/include -I$(PROTO)/include  -I. -I../src -I../../common
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
	rm -f *.o *.exe *.d ../../common/*.o

print-%: ; @echo $* = $($*)

-include $(DEPENDS)

vars:
	$(foreach v,$(filter-out $(VARS_OLD) VARS_OLD,$(.VARIABLES)), $(info $(v) = $($(v))))

