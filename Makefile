# This is a makefile.
# This makes a parallel simulation for different dielectric problem
# Use option -p in CC for profiling with gprof

PROG = np_electrostatics_lab
OBJ = main.o interface.o functions.o parallel_precal.o pfmdforces.o pcpmdforces.o penergies.o fmd.o cpmd.o

CC = g++ -O3 -g -fopenmp -Wall

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

main.o: utility.h interface.h particle.h vertex.h BIN.h control.h functions.h precalculations.h thermostat.h fmd.h 
interface.o: interface.h functions.h
functions.o: functions.h fmd.h
parallel_precal.o: precalculations.h
pfmdforces.o: forces.h
pcpmdforces.o: forces.h
penergies.o: energies.h
fmd.o: fmd.h
cpmd.o: particle.h vertex.h interface.h thermostat.h control.h forces.h energies.h functions.h

clean:
	rm -f *.o

dataclean:
	rm -f outfiles/*.dat outfiles/*.xyz  outfiles/*.lammpstrj  datafiles/*.dat verifiles/*.dat computedfiles/*.dat
