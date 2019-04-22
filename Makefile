#This make file builds the sub folder make files
PROG = np_electrostatics_lab
JOBSCR = iu_cluster_job_script.pbs
TESTDiskSCR = test_disk.pbs
TESTSphereSCR = test_sphere.pbs
BIN = bin
BASE = src
SCRIPT = scripts

all:
	@echo "Starting build of the $(BASE) directory";
ifeq ($(CCF),BigRed2)	
	+$(MAKE) -C $(BASE) cluster-install
else ifeq ($(CCF),nanoHUB)
	+$(MAKE) -C $(BASE) install
else
	+$(MAKE) -C $(BASE) local-install
endif
	@echo "Ending the build of the $(BASE) directory";
	@echo "installing the $(PROG) into $(BIN) directory"; cp -f $(BASE)/$(PROG) $(BIN)

install: all
	create-dirs

cluster-install: create-dirs
	make CCF=BigRed2 all


nanoHUB-install: create-dirs
	make CCF=nanoHUB all

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(BIN) directory"
	if ! [ -d $(BIN) ]; then mkdir $(BIN); fi
	if ! [ -d $(BIN)/outfiles ]; then mkdir $(BIN)/outfiles; fi
	if ! [ -d $(BIN)/datafiles ]; then mkdir $(BIN)/datafiles; fi
	if ! [ -d $(BIN)/verifiles ]; then mkdir $(BIN)/verifiles; fi
	if ! [ -d $(BIN)/computedfiles ]; then mkdir $(BIN)/computedfiles; fi
	@echo "Directory creation is over."

cluster-submit:
	@echo "Installing jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(JOBSCR) $(BIN)
	+$(MAKE) -C $(BIN) submit

cluster-test-disk-submit:
	@echo "Installing test jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(TESTDiskSCR) $(BIN)
	+$(MAKE) -C $(BIN) test_disk

cluster-test-sphere-submit:
	@echo "Installing test jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(TESTSphereSCR) $(BIN)
	+$(MAKE) -C $(BIN) test_sphere

clean:
	rm -f $(BASE)/*.o
	rm -f $(BASE)/$(PROG)
	rm -f $(BIN)/$(PROG)

dataclean:
	rm -f $(BIN)/outfiles/*.dat $(BIN)/outfiles/*.xyz  $(BIN)/outfiles/*.lammpstrj  $(BIN)/datafiles/*.dat verifiles/*.dat $(BIN)/computedfiles/*.dat
	rm -f $(BIN)/*.log
	rm -f $(BIN)/*.pbs

.PHONY: all clean
