#==============================================================================
#  Makefile for SanSEIR (GCC) 
#
#  Author:  Scott Collis
#
#  Revised: 3-29-2020 
#
#==============================================================================
NAME     = sanseir.exe
DEBUG    =
OPT      = -O2
FFLAGS   = -cpp -ffixed-line-length-120 -freal-4-real-8 -fdefault-real-8 \
           -fdefault-integer-8 -std=legacy $(DEFINES) $(OPT) $(DEBUG)
F90FLAGS = -cpp -freal-4-real-8 -fdefault-real-8 -fdefault-integer-8 \
           $(DEFINES) $(OPT) $(DEBUG)
OFLAGS   = $(OPT) $(DEBUG) -o $(NAME)
LIB      =
FC       = gfortran
F77      = gfortran

.SUFFIXES: .f90 

OBJECTS = sanseir.o

all: $(NAME) 

docs:
	doxygen

$(NAME): $(OBJECTS)
	$(FC) $(OFLAGS) $(OBJECTS) $(LIB)

.f90.o:
	$(FC) $(F90FLAGS) -c $*.f90 

.f.o:
	$(F77) $(FFLAGS) -c $*.f

clean:
	$(RM) *.o *.mod *.exe 
