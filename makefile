#------------------------------------------------------------------------------#
#
#  File:        makefile
#
#  Description: 
#
#------------------------------------------------------------------------------#
LIBRARY  := libgreenfish.a

SRC_DIR  := ./sources
OBJ_DIR  := ./build
INC_DIR  := ./include

DEFAULT: $(LIBRARY)

#------------------------------------------------------------------------------#
# Flags and libraries for compiling
#------------------------------------------------------------------------------#
# MPI c++ compiler
CC = $(shell which mpic++)

# OpenMPI library path
OPENMPI := $(CC:%/bin/mpic++=%)

# C compiler flags
CCFLAGS  = --std=c++11 -O3

# Include paths
INCLDIRS = -I$(OBJ_DIR) -I$(INC_DIR) -I$(OPENMPI)/include

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

$(warning Setup directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -d $(INC_DIR) || mkdir $(INC_DIR))
$(shell test -e $(LOG) && rm $(LOG))

# Dont delete the given intermediate files
.SECONDARY:

#------------------------------------------------------------------------------#
# Preprocess
#------------------------------------------------------------------------------#
CPCMD = $(CPP) $(CPPFLAGS) $(INCLDIRS) $< $@

$(OBJ_DIR)/%.cpp : $(SRC_DIR)/%.cpp
	@printf "  PRE-PROCESSING   %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#------------------------------------------------------------------------------#
# Compile
#------------------------------------------------------------------------------#
COMPILECMD = $(CC) $(CCFLAGS) $(LIBDIRS) $(INCLDIRS) -c -o $@ $< $(LIBS)

$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.cpp
	@printf "  COMPILING        %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#------------------------------------------------------------------------------#
# Link
#------------------------------------------------------------------------------#
ARCMD = ar crs $@ $(OBJECTS)

$(LIBRARY): $(OBJECTS)
	@printf "  LINKING          %-42s" "Creating library archive"; \
	$(RUN) "$(ARCMD)" $(LOG) "Error Creating Archive"; \
	$(shell cp $(SRC_DIR)/*.hpp $(INC_DIR) 2> $(LOG))

#------------------------------------------------------------------------------#
# Explicit dependencies
#------------------------------------------------------------------------------#
$(OBJ_DIR)/class_greenfish.cpp:\
  $(SRC_DIR)/class_greenfish.hpp\
  $(SRC_DIR)/class_partition.hpp\
  $(SRC_DIR)/class_communication.hpp\
  $(SRC_DIR)/greenfish/*.cpp
$(OBJ_DIR)/class_greenfish.o:\
  $(OBJ_DIR)/class_partition.o\
  $(OBJ_DIR)/class_communication.o\
  $(OBJ_DIR)/class_pencil.o
#A
#B
#C
$(OBJ_DIR)/class_communication.cpp:\
  $(SRC_DIR)/class_communication.hpp\
  $(SRC_DIR)/class_partition.hpp\
  $(SRC_DIR)/communication/*.cpp
$(OBJ_DIR)/class_communication.o:
#D
#E
#F
$(OBJ_DIR)/class_pencil.cpp:\
  $(SRC_DIR)/class_pencil.hpp\
  $(SRC_DIR)/pencil/*.cpp
$(OBJ_DIR)/class_pencil.o:
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
$(OBJ_DIR)/class_partition.cpp:\
  $(SRC_DIR)/class_partition.hpp\
  $(SRC_DIR)/partition/*.cpp
$(OBJ_DIR)/class_partition.o:
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

