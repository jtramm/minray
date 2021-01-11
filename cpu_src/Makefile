#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = no
OPENMP      = yes
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = minray

source = \
main.c \
simulation.c \
ray_trace_kernel.c \
flux_attenuation_kernel.c \
update_isotropic_sources_kernel.c \
normalize_scalar_flux_kernel.c \
add_source_to_scalar_flux_kernel.c \
compute_cell_fission_rates_kernel.c \
rand.c \
init.c \
io.c \
utils.c \
neighbor_list_a.c

header = \
minray.h \
neighbor_list_a.h

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm

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
	vim -p $(source) $(header)

run:
	./$(program)
