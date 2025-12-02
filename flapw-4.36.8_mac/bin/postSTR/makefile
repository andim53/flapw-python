include ../../MAKEFILES/make.conf
include ../postLIB/post.conf

.SUFFIXES:
.SUFFIXES: .o .f90 .f

SHELL=/bin/sh
PROGRAM=str
CWD=$(shell pwd)

flapw=pflapw
fexist = $(shell ls ../../ | grep ${flapw})
ifeq (${fexist},${flapw})
  MODI = ${MODIP}
else
  MODI = ${MODIS}
endif

MODS = mod_parameters.o mod_str.o mod_end.o

OBJS = mainstr.o read_lapwin.o read_strin.o gen_pos.o get_posdist.o\
       wrt_lapwin.o wrt_str.o ang2mat.o

.f.o:
		$(CF) $(FFLAGS) $(IPATH)../../MODI -c $(<F)

main:	$(PROGRAM)

$(PROGRAM):	$(MODS) $(OBJS)
		$(CF) $(LDFLAGS) -o x$(PROGRAM) $(MODI) $(MODS) $(OBJS) $(PLIB) $(TIME) $(FLIB) $(VERS) $(LIBS)

clean:
		-rm -f x$(PROGRAM) $(OBJS) $(MODS) *.mod *_genmod.f90
