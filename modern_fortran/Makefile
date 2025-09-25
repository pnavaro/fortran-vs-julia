SRCS =	solveur_yee.f08 sorties.f08 commun.f08
OBJS =	solveur_yee.o sorties.o commun.o
F90 = mpif90
OPT = -O3 -mtune=native -march=native
F90FLAGS = -cpp $(OPT) 
LDFLAGS = $(OPT) 

all: maxwell_serial_fortran

maxwell_serial_fortran: maxwell_serial.o 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $^ $(LIBS)

clean:
	rm -rf $(PROG) *.o *.mod data/* *.gnu core error.* output.* *.a

.SUFFIXES: $(SUFFIXES) .f08

.f08.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f08

commun.o: commun.f08
solveur_yee.o: solveur_yee.f08 commun.o
sorties.o: sorties.f08 commun.o
maxwell_serial.o: maxwell_serial.f08 $(OBJS)
