# <img src="https://github.com/WMD-group/Eris/blob/master/eris.jpg" width="48"> Eris
Codes to simulate and analyse on-lattice cation disorder in kesterite-structured Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS)

## Getting started
To run Eris locally, clone this repository and run `make` in the terminal to compile. 
To run Eris simply type `./eris`. All simulation parameters (apart from system size) are set by flags in eris.cfg. To set the system size, you need to edit the number after `#define X,Y,Z` in eris-config.c to set the X, Y and Z lattice dimensions respectively. It is necessary to recompile Eris after changing the lattice dimensions, but it is NOT necessary to recompile after changing flags in eris.cfg.

Eris can also be run in parallel over different temperatures on your local machine, for this use `make parallel`. But note that for this case, the temperature range and increments are *not* set by eris.cfg, it is now set in the Makefile by `seq 0 50 1000 | parallel ./eris {}`, which in this example corresponds to running from T=0K to T=1000K in step sizes of 50K.

## Eris on HPC
It is also possible to run Eris on a high-performance computing system. 
To run on cx1 at Imperial College, use files contained in the directory `HPC_parallel_cx1`. Files contained in `HPC_src` here are used to compile Eris on cx1. It is necessary to have libconfig in the directory here as it is linked in the Makefile. To compile Eris on cx1 type `module load intel-suite` (or `module load intel` on the Hartree Centre machines Iden and Napier), followed by `make cx1-icc`. The folder 'Eris_run_files' contains just files needed to run eris, which are: eris.cfg to set simulation parameters and a submission script. Note that for HPC eris, simular to running in parallel locally, the temperature range and increments are no longer set in eris.cfg. It is instead set in the submission script 'submit_eris.pbs' on the same line as the executable, where you must set the location of the Eris binary you have compiled on cx1 and wish to use. (Note again, that the system size must be set when compiling Eris)

## Analysing data
Bash and python scripts used to analyse data are contained in 'analysis-scripts/final_setup', where there are also some slides explaining the steps for each type of analysis. 
However, there are now ipython notebooks for the same analysis scripts with explanations and sample data (in the 'analysis-notebooks' directory.

## Publication

The underlying physical model and applications are discussed in:

* **Thermodynamically limited Cu-Zn order in Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS) from Monte Carlo simulations**  S. K. Wallace, J. M. Frost and A. Walsh (2018)

## Contributors
Jarvist Moore Frost and Suzanne K. Wallace. Codes started by Jarvist Moore Frost; September 19th 2014. Based on prior [Amphisbaena](https://github.com/jarvist/Amphisbaena) and [Starrynight](https://github.com/WMD-group/StarryNight) codes
