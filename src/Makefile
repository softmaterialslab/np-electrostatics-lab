# This is a makefile.
# This makes a parallel simulation for different dielectric problem
# Use option -p in CC for profiling with gprof

PROG = np_electrostatics_lab
OBJ = main.o interface.o functions.o parallel_precal.o pfmdforces.o pcpmdforces.o penergies.o fmd.o cpmd.o

CC = CC -O3 -g -fopenmp -Wall

LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization

CFLAG = -c

OFLAG = -o

all: $(PROG)

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)

%.o : %.cpp
	$(CC) -c $(LFLAG) $< -o $@

main.o: functions.h precalculations.h
interface.o: interface.h functions.h
functions.o: functions.h
parallel_precal.o: precalculations.h
pfmdforces.o: forces.h
pcpmdforces.o: forces.h
penergies.o: energies.h
fmd.o: functions.h
cpmd.o: functions.h
