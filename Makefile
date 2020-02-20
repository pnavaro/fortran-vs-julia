PROG = maxyee maxyee_par

SRCS =	solveur_yee.F90 sorties.F90 commun.F90 
OBJS =	solveur_yee.o sorties.o commun.o
F90 = mpif90
OPT = -O3
F90FLAGS = -cpp $(OPT) 
LDFLAGS = $(OPT) 

all: $(PROG)

maxyee: main.o 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) main.o $(LIBS)

maxyee_par: main_mpi.o 
	$(F90) $(LDFLAGS) -o $@ $(OBJS) main_mpi.o $(LIBS)

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
main.o: main.F90 $(OBJS)
main_mpi.o: main_mpi.F90 $(OBJS)
