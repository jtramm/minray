#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = minray

source = \
main.c \
simulation.c \
rand.c \
init.c \
io.c \
utils.c \
cl_utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm -lOpenCL

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3 -flto
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# OpenMP Flags
ifeq ($(OPENMP),yes)
  CFLAGS += -fopenmp -DOPENMP
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) minray.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c minray.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) *.cl minray.h

run:
	./$(program)
