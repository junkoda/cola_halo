CC   ?= gcc
EXEC = rescale_power density_grid subsampled_particles
all: $(EXEC)


GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(GSL_DIR)

LIBS := -lm -lgsl -lgslcblas

CFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS   += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)


OBJS1 = rescale_power.o 

rescale_power: $(OBJS1)
	$(CC) $(OBJS1) $(LIBS) -o $@


OBJS2 = density_grid.o
density_grid: $(OBJS2)
	$(CC) $(OBJS2) -o $@

OBJS3 = subsampled_particles.o
subsampled_particles: $(OBJS3)
	$(CC) $(OBJS3) -o $@



.PHONY: clean

clean:
	rm -f $(EXEC) $(OBJS1) $(OBJS2) $(OBJS3)
