# Makefile for comb
# Jason McEwen


# ======== OPTIONS ========

USEPGPLOT = no


# ======== COMPILER ========

FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif
OPT = $(OPTPGPLOT) -m64


# ======== LINKS ========

PROGDIR = /Users/mcewen/src

HPIXDIR = $(PROGDIR)/Healpix
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2DOC  = $(S2DIR)/doc

COMBDIR = $(PROGDIR)/comb
COMBLIB = $(COMBDIR)/lib
COMBLIBNM = comb
COMBINC = $(COMBDIR)/include
COMBSRC = $(COMBDIR)/src/mod
COMBPROG = $(COMBDIR)/src/prog
COMBBIN = $(COMBDIR)/bin
COMBDOC = $(COMBDIR)/doc

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC) -I$(COMBINC)


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

LDFLAGS =  -L$(COMBLIB) -l$(COMBLIBNM) \
	   -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) \
           -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
           $(LDFLAGSPGPLOT)


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

COMBOBJ = $(COMBINC)/comb_error_mod.o \
          $(COMBINC)/comb_obj_mod.o   \
          $(COMBINC)/comb_csky_mod.o  \
          $(COMBINC)/comb_tmplmap_mod.o \
          $(COMBINC)/comb_tmplalm_mod.o


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:     $(COMBLIB)/lib$(COMBLIBNM).a 

prog:    $(COMBBIN)/comb_csim $(COMBBIN)/comb_objgen

$(COMBINC)/%.o: $(COMBSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(COMBINC)

$(COMBINC)/%.o: $(COMBPROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 


# Library

$(COMBLIB)/lib$(COMBLIBNM).a: $(COMBOBJ)
	ar -r $(COMBLIB)/lib$(COMBLIBNM).a $(COMBOBJ)


# Documentation

docs:
	f90doc $(COMBSRC)/*.f90
	f90doc $(COMBPROG)/*.f90
	ln_multi $(S2DOC)/s2_*
	ln_multi $(S2DOC)/index_*
	mv *.html $(COMBDOC)/.
	addstyle $(COMBDOC)/comb_*

cleandocs:
	rm -f $(COMBDOC)/comb_*.html
	rm -f $(COMBDOC)/s2_*.html
	rm -f $(COMBDOC)/index_s2.html


# Clean up

clean:	tidy
	rm -f $(COMBINC)/*.mod
	rm -f $(COMBINC)/*.o
	rm -f $(COMBLIB)/lib$(COMBLIBNM).a
	rm -f $(COMBBIN)/*

tidy:
	rm -f *.mod
	rm -f $(COMBSRC)/*~
	rm -f $(COMBPROG)/*~


# Module dependencies

$(COMBINC)/comb_error_mod.o: $(COMBSRC)/comb_error_mod.f90
$(COMBINC)/comb_tmplmap_mod.o:  $(COMBSRC)/comb_tmplmap_mod.f90
$(COMBINC)/comb_tmplalm_mod.o:  $(COMBSRC)/comb_tmplalm_mod.f90
$(COMBINC)/comb_obj_mod.o:   $(COMBSRC)/comb_obj_mod.f90  \
                               $(COMBINC)/comb_error_mod.o
$(COMBINC)/comb_csky_mod.o:  $(COMBSRC)/comb_csky_mod.f90 \
                               $(COMBINC)/comb_obj_mod.o  \
                               $(COMBINC)/comb_error_mod.o


# Program dependencies and compilation

$(COMBBIN)/comb_csim.o:	$(COMBPROG)/comb_csim.f90 lib
$(COMBBIN)/comb_csim:      $(COMBINC)/comb_csim.o 
	$(FC) -o $(COMBBIN)/comb_csim $(COMBINC)/comb_csim.o \
	$(LDFLAGS) $(PPFLAGS)

$(COMBBIN)/comb_objgen.o:	$(COMBPROG)/comb_objgen.f90 lib
$(COMBBIN)/comb_objgen:      $(COMBINC)/comb_objgen.o 
	$(FC) -o $(COMBBIN)/comb_objgen $(COMBINC)/comb_objgen.o \
	$(LDFLAGS) $(PPFLAGS)


