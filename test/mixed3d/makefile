#------------------------------------------------------------------------------#
#
#  File:        makefile
#
#  Description: 
#
#------------------------------------------------------------------------------#
PROGRAM  := mixed3d

MAINFILE := ./main.cpp

SRC_DIR  := ./sources
OBJ_DIR  := ./build

default: $(PROGRAM)

#------------------------------------------------------------------------------#
# Dynamic path to Poisson solver library
#------------------------------------------------------------------------------#
POISSON = ./../../

#------------------------------------------------------------------------------#
# Flags and libraries for compiling
#------------------------------------------------------------------------------#
# MPI c++ compiler
CC = $(shell which mpic++)

# OpenMPI library path
OPENMPI = $(CC:%/bin/mpic++=%)

# c++ compiler flags
CCFLAGS   = --std=c++11 -O3

# Library flags
LIBDIRS = -L$(POISSON)

LIBS = -lpoissonsolver

# Include paths
INCLDIRS = -I$(POISSON)/include -I$(OPENMPI)/include

# Preprocessor compiler
CPP = cpp

# Preprocessor compiler flags
CPPFLAGS = --std=c++11 -O3 -D__verb

#------------------------------------------------------------------------------#
# write scripts
#------------------------------------------------------------------------------#
SHELL := /bin/sh
LOG   := compile.log
RUN   := $(SHELL) ./etc/make_output.sh

#------------------------------------------------------------------------------#
# Setup build directories
#------------------------------------------------------------------------------#
NAMES   := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
SOURCES := $(NAMES:%.cpp=$(SRC_DIR)/%.cpp) 
OBJECTS := $(NAMES:%.cpp=$(OBJ_DIR)/%.o) 
#$(OBJ_DIR)/main.o

$(warning Setup directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -d output || mkdir output)
$(shell test -e $(LOG) && rm $(LOG))


# Dont delete the given intermediate files
.SECONDARY:

#------------------------------------------------------------------------------#
# Preprocess
#------------------------------------------------------------------------------#
CPCMD = $(CPP) $(CPPFLAGS) $(INCLDIRS) $< $@
$(OBJ_DIR)/%.cpp : $(SRC_DIR)/%.cpp
	@printf "  CPP  %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#$(OBJ_DIR)/main.cpp : $(MAINFILE)
#	@printf "  CPP  %-42s" $<; \
#	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#------------------------------------------------------------------------------#
# Compile
#------------------------------------------------------------------------------#
COMPILECMD = $(CC) $(CCFLAGS) $(INCLDIRS) $(LIBDIRS) -c -o $@ $< $(LIBS)
$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.cpp
	@printf "  CC   %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#------------------------------------------------------------------------------#
# Link
#------------------------------------------------------------------------------#
LINKCMD = $(CC) $(CCFLAGS) $(LIBDIRS) $(INCLDIRS) -o $@ $(OBJECTS) $(LIBS)
$(PROGRAM): $(OBJECTS)
	@printf " LINK  %-42s" $@; \
	$(RUN) "$(LINKCMD)" $(LOG) "Linking Error"

#------------------------------------------------------------------------------#
# Explicit dependencies
#------------------------------------------------------------------------------#
$(OBJ_DIR)/main.cpp:
# $(SRC_DIR)/*.hpp
$(OBJ_DIR)/main.o:\
  $(POISSON)/*.a
#A
#B
#C
#D
#E
#F
#G
#H
#I
#J
#K
#L
#M
#N
#O
#P
#Q
#R
#S
#T
#U
#V
#W
#X
#Y
#Z

#------------------------------------------------------------------------------#
# Clean
#------------------------------------------------------------------------------#
clean:
	rm -fr $(OBJ_DIR)
	rm -f  $(PROGRAM)
	rm -f  $(LOG)

run: $(PROGRAM)
	mpirun -np 4 $(PROGRAM)
