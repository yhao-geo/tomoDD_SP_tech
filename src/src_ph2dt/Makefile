CMD	= ph2dtN3
CC	= gcc
FC	= gfortran

SRCS	= $(CMD).f \
	  cluster.f datetime.f delaz.f ifindi.f \
	  indexx.f indexxi.f sorti.f trimlen.f
CSRCS	= atoangle_.c atoangle.c rpad_.c sscanf3_.c
OBJS	= $(SRCS:%.f=%.o) $(CSRCS:%.c=%.o)
INCLDIR	= ./include

CFLAGS	= -I$(INCLDIR) -O -mcmodel=large
FFLAGS	= -O -I$(INCLDIR) -std=legacy -fno-range-check -Wno-argument-mismatch -mcmodel=large

all: $(CMD)

$(CMD): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(CMD) *.o core a.out junk

# Include file dependencies
atoangle_.o : $(INCLDIR)/compat.h
atoangle_.o : $(INCLDIR)/f77types.h
ph2dtN3.o   : $(INCLDIR)/ph2dt.inc
