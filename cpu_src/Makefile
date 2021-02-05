#===============================================================================
# User Options
#===============================================================================

COMPILER       = gnu
OPTIMIZE       = yes
DEBUG          = no
OPENMP         = yes
PROFILE        = no
ifeq ($(NL_ALGORITHM),)
  NL_ALGORITHM = c
endif
QUIET          = yes

#===============================================================================
# Sets Basic Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm

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
utils.c 

header = \
minray.h \
atomic_wrapper.h

ifeq ($(NL_ALGORITHM),a)
  source += neighbor_list_a.c
	header += neighbor_list_a.h
	CFLAGS += -DALGORITHM_A
endif
ifeq ($(NL_ALGORITHM),b)
  source += neighbor_list_b.c
	header += neighbor_list_b.h
	CFLAGS += -DALGORITHM_B
endif
ifeq ($(NL_ALGORITHM),c)
  source += neighbor_list_c.c
	header += neighbor_list_c.h
	CFLAGS += -DALGORITHM_C
endif
ifeq ($(NL_ALGORITHM),d)
  source += neighbor_list_d.c
	header += neighbor_list_d.h
	CFLAGS += -DALGORITHM_D
endif
ifeq ($(NL_ALGORITHM),e)
  source += neighbor_list_e.c
	header += neighbor_list_e.h
	CFLAGS += -DALGORITHM_E
endif
ifeq ($(NL_ALGORITHM),f)
  source += neighbor_list_f.c
	header += neighbor_list_f.h
	CFLAGS += -DALGORITHM_F
endif
ifeq ($(NL_ALGORITHM),g)
  source += neighbor_list_g.c
	header += neighbor_list_g.h
	CFLAGS += -DALGORITHM_G
endif
ifeq ($(NL_ALGORITHM),h)
  source += neighbor_list_h.c
	header += neighbor_list_h.h
	CFLAGS += -DALGORITHM_H
endif
ifeq ($(NL_ALGORITHM),i)
  source += neighbor_list_i.c
	header += neighbor_list_i.h
	CFLAGS += -DALGORITHM_I
endif
ifeq ($(NL_ALGORITHM),j)
  source += neighbor_list_j.c
	header += neighbor_list_j.h
	CFLAGS += -DALGORITHM_J
endif
ifeq ($(NL_ALGORITHM),k)
  source += neighbor_list_k.c
	header += neighbor_list_k.h
	CFLAGS += -DALGORITHM_K -DLEVEL_0
endif
ifeq ($(NL_ALGORITHM),l)
  source += neighbor_list_k.c
	header += neighbor_list_k.h
	CFLAGS += -DALGORITHM_K -DLEVEL_3
endif

obj = $(source:.c=.o)

#===============================================================================
# Misc Flags
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

# Quiet the iteration status
ifeq ($(QUIET),yes)
  CFLAGS += -DQUIET
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
