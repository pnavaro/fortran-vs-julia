SRCS =	solveur_yee.F90 sorties.F90 commun.F90 
OBJS =	solveur_yee.o sorties.o commun.o
F90 = mpif90
OPT = -O0 -g -fbounds-check -O3
F90FLAGS = -cpp $(OPT) 
LDFLAGS = $(OPT) 

all: maxwell_serial_fortran maxwell_mpi_fortran

maxwell_serial_fortran: maxyee.o 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $^ $(LIBS)

maxwell_mpi_fortran: maxyee_mpi.o 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $^ $(LIBS)

clean:
	rm -rf $(PROG) *.o *.mod data/* *.gnu core error.* output.* *.a

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.F90

commun.o: commun.F90
solveur_yee.o: solveur_yee.F90 commun.o
sorties.o: sorties.F90 commun.o
maxyee.o: maxyee.F90 $(OBJS)
maxyee_mpi.o: maxyee_mpi.F90 $(OBJS)
