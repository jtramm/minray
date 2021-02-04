#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
ALGORITHM   = b

#===============================================================================
# Program name & source code list
#===============================================================================

program = minray

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm -lOpenCL

source = \
main.c \
simulation.c \
rand.c \
init.c \
io.c \
utils.c \
cl_utils.c

header = \
minray.h \
parameters.h \
atomic_float.h

obj = $(source:.c=.o)

ifeq ($(ALGORITHM),b)
  source += neighbor_list_host_b.c
	header += neighbor_list_b.h
	CFLAGS += -DALGORITHM_B
endif

ifeq ($(ALGORITHM),h)
  source += neighbor_list_host_h.c
	header += neighbor_list_h.h
	CFLAGS += -DALGORITHM_H
endif

ifeq ($(ALGORITHM),j)
  source += neighbor_list_host_j.c
	header += neighbor_list_j.h
	CFLAGS += -DALGORITHM_J
endif

ifeq ($(ALGORITHM),k)
  source += neighbor_list_host_k.c
	header += neighbor_list_k.h
	CFLAGS += -DALGORITHM_K
endif

#===============================================================================
# Sets Flags
#===============================================================================


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

$(program): $(obj) $(header) Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c $(header) Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) *.cl $(header)

run:
	./$(program)
