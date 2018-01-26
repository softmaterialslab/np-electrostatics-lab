# This is a makefile.
# This makes a parallel simulation for different dielectric problem
# Use option -p in CC for profiling with gprof

PROG = np_electrostatics_lab
OBJ = main.o interface.o functions.o parallel_precal.o pfmdforces.o pcpmdforces.o penergies.o fmd.o cpmd.o

CC = g++ -O3 -g -fopenmp -Wall -g

LFLAG = -lgsl -lgslcblas -lm -lboost_program_options

CFLAG = -c

OFLAG = -o

all: $(PROG)

#dirs: 
#	@echo "Creating needed sub-directories in the current directory"; mkdir outfiles; mkdir datafiles; mkdir verifiles; mkdir computedfiles

install: all
	@echo "Installing $(PROG) into current directory"

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

clean:
	rm -f *.o

dataclean:
	rm -f outfiles/*.dat outfiles/*.xyz  outfiles/*.lammpstrj  datafiles/*.dat verifiles/*.dat computedfiles/*.dat
