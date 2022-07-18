# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += cgebhart

## Define the variables needed by Make.example

PROTO_HOME = $(HOME)/gitroot/proto
# PROTO_HOME = $(SCRATCH)/proto
# PROTO_HOME = $(HOME)/proto

# the base name(s) of the application(s) in this directory
ebase = test imexTest trigTest dipoleTest waveTest testMPI testMove

# the location of the Chombo "lib" directory
CHOMBO_HOME =$(HOME)/svnroot/Chombo/lib
# CHOMBO_HOME =$(SCRATCH)/Chombo/lib
# CHOMBO_HOME =$(HOME)/Chombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTools BoxTools BaseTools
# the locations of the source code directories
base_dir = .
src_dirs = 

# input file for 'run' target
INPUT =
# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# custom variables

# application-specific variables
CPPFLAGS += -DUSE_PROTO=TRUE -DDIM=$(DIM) -I./../src -I$(PROTO_HOME)/include
# application-specific targets

print-% : ; @echo $* = $($*)

dataclean : ; rm *.hdf5 *.err *.out
