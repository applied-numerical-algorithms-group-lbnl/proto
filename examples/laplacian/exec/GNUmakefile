VARS_OLD := $(.VARIABLES)
DEPXX=g++
CXX=g++
DIM=2
PROTO=../../..

CXXOPTFLAGS=-g 
TARGET:=laplacian.exe

SRC=$(patsubst ../src/%,%,$(wildcard ../src/*.cpp)) $(wildcard *.cpp) $(wildcard ../../common/*.cpp) 
NODIRSRC = $(notdir $(SRC))
NODIROBJ = $(subst .cpp,.o,$(NODIRSRC))
OBJ=$(patsubst %.o,o/%.o, $(NODIROBJ))

NODIRDEP = $(subst .cpp,.d,$(NODIRSRC))
DEPENDS=$(patsubst %.d,d/%.d, $(NODIRDEP))

VPATH=o. ../src  ../../common


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



