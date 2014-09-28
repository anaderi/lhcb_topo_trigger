LD      = g++
CC      = gcc
FLAGS   = -O2 -Wall
INCLUDE = -I /usr/local/include/root
# ROOT_FLAGS      = $(shell $(ROOTSYS)/bin/root-config --cflags)
# ROOT_LIBS       = $(shell $(ROOTSYS)/bin/root-config --libs)
LIBS = -L /usr/local/lib/root -lTree -lPhysics
# -lRooFit -lRooStats
# FLAGS += ${ROOT_FLAGS}
#
# general pattern to build executables from object files
#

all: generate_csv

%: %.cc
	$(LD) $(FLAGS) $(INCLUDE) -c -o $*.o $*.cc
	$(LD) $(FLAGS) $*.o $(LIBS) -o $*
