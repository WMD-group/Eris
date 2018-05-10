# <img src="https://github.com/WMD-group/Eris/blob/master/eris.jpg" width="48"> Eris
Monte Carlo code to simulate and analyse on-lattice cation disorder in kesterite-structured Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS)

## Getting started
To run Eris locally, clone this repository and run `make` in the terminal to compile. 
To run Eris simply type `./eris`. All simulation parameters are set by flags in eris.cfg and details of flag settings for certain analyses are included in the corresponding ipython notebook in the 'analysis-notebooks' directory.

Eris can also be run in parallel over different temperatures on your local machine, for this use `make parallel`. But note that for this case, the temperature range and increments are **not** set by eris.cfg, it is now set in the Makefile by `seq 0 50 1000 | parallel ./eris {}`, which in this example corresponds to running from T=0K to T=1000K in step sizes of 50K.

See the wiki page for more information.

## Eris on HPC
It is also possible to run Eris on a high-performance computing system. 
To run on cx1 at Imperial College, use `make cx1-icc` to compile. It is necessary to have libconfig as it is linked in the Makefile. To compile Eris on cx1 type `module load intel-suite` (or `module load intel` on the Hartree Centre machines Iden and Napier), followed by `make cx1-icc`. The following executable line is required in the submission script for the HPC system: `seq 0 50 1000 | parallel -j NUM-CORES PATH-TO-BINARY/eris {} > eris-parallel.dat`, where 'NUM-CORES' and 'PATH-TO-BINARY' need to be modified.

## Analysing data
ipython notebooks used to analyse data are contained in 'analysis-notebooks' and data is uploaded to the zenodo repository, doi: ????

## Publication

The underlying physical model and applications are discussed in:

* **Thermodynamically limited Cu-Zn order in Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS) from Monte Carlo simulations**  S. K. Wallace, J. M. Frost and A. Walsh (2018)

## Contributors
Jarvist Moore Frost and Suzanne K. Wallace. Codes started by Jarvist Moore Frost; September 19th 2014. Based on prior [Amphisbaena](https://github.com/jarvist/Amphisbaena) and [Starrynight](https://github.com/WMD-group/StarryNight) codes
