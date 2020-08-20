#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = yes
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = minray

source = \
main.c \
ray_trace_kernel.c \
rand.c \
init.c \
io.c \
clutils.c

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
  CFLAGS += -lOpenCL
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
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
	vim -p $(source) ray_trace_kernel.cl flux_kernel.cl minray.h

run:
	./$(program)
