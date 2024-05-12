# Makefile for compiling PETSc program

# Define the compiler
CC=gcc

# Define the executable name
EXE=hello

# Source files
SRC=hello.c

# PETSc installation directory and architecture
PETSC_DIR=/home/ber0061/software/petsc
PETSC_ARCH=arch-linux-c-opt

# Include and library paths
INCLUDES=-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIBS=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

# Compiler flags
CFLAGS=$(INCLUDES)
LFLAGS=$(LIBS)

# Build rule for the executable
all: $(EXE)

$(EXE): $(SRC)
	$(CC) -o $@ $< $(CFLAGS) $(LFLAGS)

# Rule to run the program with environment variable set
run: $(EXE)
	LD_LIBRARY_PATH=$(PETSC_DIR)/$(PETSC_ARCH)/lib:$$LD_LIBRARY_PATH ./$(EXE)

# Clean rule to remove compiled objects and executable
clean:
	rm -f $(EXE)

.PHONY: all run clean
