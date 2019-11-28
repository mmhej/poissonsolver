#------------------------------------------------------------------------------#
#
#  File:        makefile
#
#  Description: Make file for compiling the Poisson solver library
#
#------------------------------------------------------------------------------#
LIBRARY  := libpoissonsolver.a

SRC_DIR  := ./source
OBJ_DIR  := ./build
INC_DIR  := ./include

DEFAULT: $(LIBRARY)

#------------------------------------------------------------------------------#
# Flags and libraries for compiling
#------------------------------------------------------------------------------#
# MPI Fortran compiler
FC = $(shell which mpif90)

# OpenMPI library path
OPENMPI := $(FC:%/bin/mpif90=%)

# FC compiler flags
FCFLAGS  = -O3 -ffree-form -ffree-line-length-none
#FCFLAGS  = -O3 -ffree-form -ffree-line-length-none -fbounds-check

# Library flags
LIBDIRS = 

LIBS = 

# Include paths
INCLDIRS = -I$(INC_DIR)

# Preprocessor compiler
CPP = cpp

# Preprocessor compiler flags
CPPFLAGS = -D__double

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
$(shell test -d $(INC_DIR) || mkdir $(INC_DIR))
$(shell test -e $(LOG) && rm $(LOG))

# Dont delete the given intermediate files
.SECONDARY:

#------------------------------------------------------------------------------#
# Preprocess
#------------------------------------------------------------------------------#
CPCMD = $(CPP) $(CPPFLAGS) $< $@
$(OBJ_DIR)/%.f : $(SRC_DIR)/%.f
	@printf "  PRE-PROCESSING   %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#------------------------------------------------------------------------------#
# Compile
#------------------------------------------------------------------------------#
COMPILECMD = $(FC) $(FCFLAGS) $(LIBDIRS) $(LIBS) $(INCLDIRS) -c -o $@ $<
$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.f
	@printf "  COMPILING        %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#------------------------------------------------------------------------------#
# Link
#------------------------------------------------------------------------------#
ARCMD = ar crs $@ $(OBJECTS)
$(LIBRARY): $(OBJECTS)
	@printf "  LINKING          %-42s" "Creating library archive"; \
	$(RUN) "$(ARCMD)" $(LOG) "Error Creating Archive"; \
	$(shell mv *.mod $(INC_DIR) >> $(LOG))

#------------------------------------------------------------------------------#
# Explicit dependencies
#------------------------------------------------------------------------------#
$(OBJ_DIR)/poisson_solver.f:\
  $(SRC_DIR)/poisson_solver/*.f\
  $(SRC_DIR)/precision.h
$(OBJ_DIR)/poisson_solver.o:\
  $(OBJ_DIR)/partition.o\
  $(OBJ_DIR)/communication.o\
  $(OBJ_DIR)/pencil.o
#A
#B
#C
$(OBJ_DIR)/communication.f:\
  $(SRC_DIR)/communication/*.f\
  $(SRC_DIR)/precision.h
$(OBJ_DIR)/communication.o:
#D
#E
#F
$(OBJ_DIR)/pencil.f:\
  $(SRC_DIR)/pencil/*.f\
  $(SRC_DIR)/precision.h
$(OBJ_DIR)/pencil.o:
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
$(OBJ_DIR)/partition.f:\
  $(SRC_DIR)/partition/*.f\
  $(SRC_DIR)/precision.h
$(OBJ_DIR)/partition.o:
#U
#V
#W
#X
#Y
#Z

#------------------------------------------------------------------------------#
# Clean
#------------------------------------------------------------------------------#
#clean:
#	rm -fr $(OBJ_DIR)
#	rm -f  $(LOG)

clean:
	rm -fr $(OBJ_DIR)
	rm -f  $(LIBRARY)
	rm -f  $(LOG)

clean_all:
	rm -fr $(OBJ_DIR)
	rm -fr $(INC_DIR)
	rm -f  $(LIBRARY)
	rm -f  $(LOG)
	cd test/unbounded2d/ && make clean
	cd test/unbounded3d/ && make clean
	cd test/periodic2d/ && make clean
	cd test/periodic3d/ && make clean
	cd test/mixed2d/ && make clean
	cd test/mixed_ffp3d/ && make clean
	cd test/mixed_ppf3d/ && make clean

test_all:
	make
	cd test/unbounded2d/ && make run
	cd test/unbounded3d/ && make run
	cd test/periodic2d/ && make run
	cd test/periodic3d/ && make run
	cd test/mixed2d/ && make run
	cd test/mixed_ffp3d/ && make run
	cd test/mixed_ppf3d/ && make run

test_one:
	make
	cd test/unbounded2d/ && make run

