# np-electrostatics-lab

## Install instructions on BigRed2
* First, git clone the project
```git clone https://github.com/softmaterialslab/np-electrostatics-lab.git```
* Then, load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Next, go to the root directory:
 ```cd np-electrostatics-lab```
* Then, install the project:
```make install```
* Fianlly, submit the job:
```make cluster-submit```
* All outputs from the simulation will be stored in the data folder when the simulation is completed.


## Necessary Modules for local run
* Load the necessary modules: module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0
* Make sure to export BOOST_LIBDIR environment variable with location to the lib directory: 
```export BOOST_LIBDIR=/opt/boost/gnu/openmpi_ib/lib/```
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU: 
```export OMP_NUM_THREADS=16```

## Local Install instructions

* First, git clone the project
```git clone https://github.com/softmaterialslab/np-electrostatics-lab.git```
* Go to np-electrostatics-lab directory:
 ```cd np-electrostatics-lab```
* You should provide the following make command to make the project. This will create the executable and Install the executable (np_electrostatics_lab) into bin directory (That is np-electrostatics-lab/bin)
 ```make local-install ```
* Next, go to the bin directory: 
 ```cd bin ```
* Now you are ready to run the executable with aprun command using the following method:
 * Spehre:
```time mpirun -np 2 -N 16 ./np_electrostatics_lab -a 2.6775 -b 14.28 -e 2 -E 78.5 -V -60 -v 1 -g 1082 -m 6 -t 0.001 -s 10000 -p 100 -f 10 -M 6 -T 0.001 -k 0.0025 -q 0.001 -L 5 -l 5 -S 10000000 -P 100000 -F 100 -X 10000 -U 1000 -Y 500000 -W 1000000 -B 0.025```
 * Disk:
 ```time mpirun -np 2 -N 16 ./np_electrostatics_lab -a 2.6775 -b 14.28 -e 2 -E 78.5 -V -60 -v 1 -g 1082 -m 6 -t 0.001 -s 10000 -p 100 -f 10 -M 6 -T 0.001 -k 0.0025 -q 0.001 -L 5 -l 5 -S 10000000 -P 100000 -F 100 -X 10000 -U 1000 -Y 500000 -W 1000000 -R 0.1 -B 0.4 -G "Disk"```


## NanoHUB app page:
* https://nanohub.org/tools/nselectrostatic/



