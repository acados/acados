# This Makefile helps to compile the acados example
# test_ocp_qp_hpmpc on microcontrollers

# Makefile to help compiling acados code in a directory
# It scans the SRCDIR to find all source files,
# then makes the object for every source file,
# after that it links all the objects to compile the project
# Note: PROJECT is the name of the main file to execute
# Author: Dang Doan
# Date: 2017.03.15

# BLASFEO-specific options
LA = REFERENCE
TARGET = GENERIC

# Compiler
CC = gcc
CFLAGS = -std=c99 -Wextra -Wno-error -O2 -fPIC -march=native
# -Wall

# Project name (binary)
PROJECT = test_ocp_qp_hpmpc

# Exclude particular source files
IGNORESRC = $(SRCDIR)/$(PROJECT).$(SRCEXT)
# main.c in the code of ESP32
IGNORESRC += $(SRCDIR)/main.c

# Directories
BUILDDIR = obj
SRCDIR = .
TARGETDIR = bin
INCDIR = include
# if all .h file are put in the same directory, use:
# INCDIR = .

# Flags
IFLAGS = -I$(SRCDIR) -I$(INCDIR)
LDFLAGS = -lm
SRCEXT = c
OBJEXT = o
BINEXT = exe

# Scan for sources and objects to compile
ALL_SRCS = $(wildcard $(SRCDIR)/*.$(SRCEXT))

SRCS = $(filter-out $(IGNORESRC), $(ALL_SRCS))
# SRCS = $(filter-out $(SRCDIR)/$(PROJECT).$(SRCEXT), $(ALL_SRCS))
OBJS = $(patsubst $(SRCDIR)/%.$(SRCEXT),$(BUILDDIR)/%.$(OBJEXT),$(SRCS))

# Target file to make
ifdef BINEXT
	TARGET = $(PROJECT).$(BINEXT)
else
	TARGET = $(PROJECT)
endif

# Libraries
# LIBS =

# Targets
all: resources $(TARGET)

resources: directories

#Make the Directories if not exist
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

# Link
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(IFLAGS) -o $(TARGETDIR)/$(TARGET) $(PROJECT).$(SRCEXT) $^ $(LDFLAGS)

# Compile, with debug option for GDB and valgrind checking
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	$(CC) $(CFLAGS) $(IFLAGS) $< ${STATIC_LIBRARIES} -c  -g -O0 -o $@

clean:
	rm $(BUILDDIR)/*.$(OBJEXT) $(TARGETDIR)/$(TARGET)

.PHONY: all resources clean
