
FFLAGS = -traceback -O3 
 FFLAGS = -O3 
#flags for gfortran
CC = mpif90

OBJECTS = constants.o parameters.o i_o.o distribsinit.o initialise.o quasilinear_terms.o \
  calc_indices.o electron_source.o nonlinsK.o emprop.o applybounds.o 

HELPERMODS = constants.o parameters.o

TARGET = code


code : $(OBJECTS) beam.o
	$(CC) $(FFLAGS) -o $(TARGET) $(OBJECTS) beam.o

constants.o:constants.f90

parameters.o : constants.o
parameters : constants

i_o.o : $(HELPERMODS)
i_o: $(HELPERMODS)

distribsinit.o:$(HELPERMODS)
distribsinit:$(HELPERMODS)

initialise.o:$(HELPERMODS)
initialise:$(HELPERMODS)

quasilinear_terms.o:$(HELPERMODS)
quasilinear_terms:$(HELPERMODS)

calc_indices.o:$(HELPERMODS)
calc_indices:$(HELPERMODS)

electron_source.o:$(HELPERMODS)
electron_source:$(HELPERMODS)

nonlinsK.o:$(HELPERMODS)
nonlinsK:$(HELPERMODS)

emprop.o:$(HELPERMODS)
emprop:$(HELPERMODS)

applybounds.o:$(HELPERMODS)
applybounds:$(HELPERMODS)

beam.o:$(OBJECTS)
beam:$(OBJECTS)

# ======================================================================
# And now the general rules, these should not require modification
#FROM Makefile-Fortran
#A Makefile for a simple Fortran project
#Copyright (C) 2009 Davide Cesari, dcesari <at> arpa.emr.it
#The sources are distributed according to the GNU GPL license.
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(CC) $(FFLAGS) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(CC) $(FFLAGS) -c $<

%.o: %.F90
	$(CC) $(FFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD ./code

veryclean: clean
	rm -f *~ $(OBJECTS)


