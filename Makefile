#This make file builds the sub folder make files
PROG = np_electrostatics_lab
JOBSCR = iu_cluster_job_script.pbs
TESTSCR = test.pbs
BIN = bin
BASE = src
DATA = data
SCRIPT = scripts

all:
	@echo "Starting build of the $(BASE) directory";
ifeq ($(CCF),BigRed2)	
	+$(MAKE) -C $(BASE) cluster-install
else ifeq ($(CCF),nanoHUB)
	+$(MAKE) -C $(BASE)
else
	+$(MAKE) -C $(BASE)
endif
	@echo "Ending the build of the $(BASE) directory";
	@echo "installing the $(PROG) into $(BIN) directory"; cp -f $(BASE)/$(PROG) $(BIN)

install: all
	create-dirs
	@echo "Installing $(PROG) into $(DATA) directory"
	cp -f $(BIN)/$(PROG) $(DATA)

cluster-install: create-dirs
	make CCF=BigRed2 all
	@echo "Installing $(PROG) into $(DATA) directory"
	cp -f $(BIN)/$(PROG) $(DATA)

nanoHUB-install: create-dirs
	make CCF=nanoHUB all
	@echo "Installing $(PROG) into $(DATA) directory"
	cp -f $(BIN)/$(PROG) $(DATA)

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(DATA) directory"
	if ! [ -d $(DATA) ]; then mkdir $(DATA); fi
	if ! [ -d $(DATA)/outfiles ]; then mkdir $(DATA)/outfiles; fi
	if ! [ -d $(DATA)/datafiles ]; then mkdir $(DATA)/datafiles; fi
	if ! [ -d $(DATA)/verifiles ]; then mkdir $(DATA)/verifiles; fi
	if ! [ -d $(DATA)/computedfiles ]; then mkdir $(DATA)/computedfiles; fi
	@echo "Directory creation is over."

cluster-submit:
	@echo "Installing jobscript into $(DATA) directory"
	cp -f $(SCRIPT)/$(JOBSCR) $(DATA)
	+$(MAKE) -C $(DATA) submit

cluster-test-submit:
	@echo "Installing test jobscript into $(DATA) directory"
	cp -f $(SCRIPT)/$(TESTSCR) $(DATA)
	+$(MAKE) -C $(DATA) test

clean:
	rm -f $(BASE)/*.o
	rm -f $(BASE)/$(PROG)
	rm -f $(BIN)/$(PROG)

dataclean:
	rm -f $(DATA)/$(PROG)
	rm -f $(DATA)/$(JOBSCR)
	rm -f $(DATA)/outfiles/*.dat $(DATA)/outfiles/*.xyz  $(DATA)/outfiles/*.lammpstrj  $(DATA)/datafiles/*.dat verifiles/*.dat $(DATA)/computedfiles/*.dat
	rm -f $(DATA)/*.log
	rm -f $(DATA)/*.pbs

.PHONY: all clean
