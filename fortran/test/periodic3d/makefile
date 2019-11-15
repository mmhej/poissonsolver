#------------------------------------------------------------------------------#
#
#  File:        makefile
#
#  Description: 
#
#------------------------------------------------------------------------------#
PROGRAM  := periodic3d

MAINFILE := ./main.f

SRC_DIR  := ./source
OBJ_DIR  := ./build

default: $(PROGRAM)

#------------------------------------------------------------------------------#
# Dynamic path to Poisson solver library
#------------------------------------------------------------------------------#
POISSON = ./../..

#------------------------------------------------------------------------------#
# Flags and libraries for compiling
#------------------------------------------------------------------------------#
# MPI Fortran compiler
FC = $(shell which mpif90)

# OpenMPI library path
OPENMPI := $(FC:%/bin/mpif90=%)

# FC compiler flags
FCFLAGS  = -O3 -ffree-form -ffree-line-length-none

# Library flags
LIBDIRS = -L$(POISSON)

LIBS = -lpoissonsolver

# Include paths
INCLDIRS = -I$(POISSON)/include -I$(OPENMPI)/include

# Preprocessor compiler
CPP = cpp

# Preprocessor compiler flags
CPPFLAGS = -D__verb

#------------------------------------------------------------------------------#
# write scripts
#------------------------------------------------------------------------------#
SHELL := /bin/sh
LOG   := compile.log
RUN   := $(SHELL) ./etc/make_output.sh

#------------------------------------------------------------------------------#
# Setup build directories
#------------------------------------------------------------------------------#
NAMES   := $(notdir $(wildcard $(SRC_DIR)/*.f))
SOURCES := $(NAMES:%.f=$(SRC_DIR)/%.f) 
OBJECTS := $(NAMES:%.f=$(OBJ_DIR)/%.o) 

$(warning Setup directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -d output || mkdir output)
$(shell test -e $(LOG) && rm $(LOG))


# Dont delete the given intermediate files
.SECONDARY:

#------------------------------------------------------------------------------#
# Preprocess
#------------------------------------------------------------------------------#
CPCMD = $(CPP) $(CPPFLAGS) $< $@
$(OBJ_DIR)/%.f : $(SRC_DIR)/%.f
	@printf "  CPP  %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#------------------------------------------------------------------------------#
# Compile
#------------------------------------------------------------------------------#
COMPILECMD = $(FC) $(FCFLAGS) $(LIBDIRS) $(LIBS) $(INCLDIRS) -c -o $@ $<
$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.f
	@printf "  CC   %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#------------------------------------------------------------------------------#
# Link
#------------------------------------------------------------------------------#
LINKCMD = $(FC) $(FCFLAGS) $(OBJECTS) $(LIBDIRS) $(LIBS) $(INCLDIRS) -o $@
$(PROGRAM): $(OBJECTS)
	@printf " LINK  %-42s" $@; \
	$(RUN) "$(LINKCMD)" $(LOG) "Linking Error"

#------------------------------------------------------------------------------#
# Explicit dependencies
#------------------------------------------------------------------------------#
$(OBJ_DIR)/main.f:\
  $(POISSON)/*.a
$(OBJ_DIR)/main.o:
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
