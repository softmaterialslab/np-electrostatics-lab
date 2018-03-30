# np-electrostatics-lab

## Install instructions on BigRed2
* git clone the project
```git clone https://github.com/softmaterialslab/np-electrostatics-lab.git```
* load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Go to root directory:
 ```cd np-electrostatic-lab```
* Next install the project:
```make install```
* Fianlly submit the job:
```make cluster-submit```
* All outputs from the simulation will be stored in the data folder when the simulation is completed.

