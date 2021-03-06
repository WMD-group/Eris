[![DOI](https://zenodo.org/badge/116799125.svg)](https://zenodo.org/badge/latestdoi/116799125)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build
Status](https://travis-ci.org/WMD-group/Eris.svg?branch=master)](https://travis-ci.org/WMD-group/Eris)

# <img src="eris.jpg" width="48"> Eris
Monte Carlo codes to simulate and analyse on-lattice cation disorder in kesterite-structured Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS)

## Getting started

Eris uses [libconfig](https://github.com/hyperrealm/libconfig) for reading the configuration file. 
On a Debian-based Linux, this can be installed by `sudo apt-get install libconfig-dev`.

To run locally, clone this repository and run `make` to compile. 
To run the compiled binary, type `./eris`. All simulation parameters are set by flags in `eris.cfg` and details of flag settings for certain analyses are included in the corresponding ipython notebook in the [analysis-notebooks](/analysis-notebooks/) directory.

Eris can be (trivially) parallelised over temperature with GNU `parallel`. Use `make parallel` for an example. Note that this overrides the temperature range and increments set by `eris.cfg`, via the command line interface. The `Makefile` example of `seq 0 50 1000 | parallel ./eris {}` corresponds to running from T=0 K to T=1000 K in steps of 50 K.

See the [Wiki](https://github.com/WMD-group/Eris/wiki/) for more information.

## Eris on HPC
Eris can be run in a queue on high-performance computing systems. 

To compile for [CX1](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/hpc/hpc-service-support/service/) at Imperial College with the intel compiler, use `make cx1-icc`. This relies on first compiling a local copy of libconfig, to statically link into the Eris binary. A similar process can be used on other machines to which you don't have administrator rights. 
This would first require the loading of the necessary Intel development environment module (CX1: `module load intel-suite`, [Hartree Centre](http://community.hartree.stfc.ac.uk) Iden/Napier: `module load intel`). The following executable line is required in the job submission script: `seq 0 50 1000 | parallel -j NUM-CORES PATH-TO-BINARY/eris {} > eris-parallel.dat`, where `NUM-CORES` and `PATH-TO-BINARY` need to be modified. Note again here that the temperature range set in `eris.cfg` will be overwritten.

## Analysing data
ipython notebooks used to analyse data are contained in [analysis-notebooks](/analysis-notebooks/). Raw data from the publications listed below has been uploaded to a [Zenodo repository](https://doi.org/10.5281/zenodo.1251122).

## Publication

The underlying physical model and applications are discussed in:

* **Atomistic insights into the order–disorder transition in Cu<sub>2</sub>ZnSnS<sub>4</sub> solar cells from Monte Carlo simulations**  S. K. Wallace, J. M. Frost and A. Walsh, DOI: 10.1039/C8TA04812F (Paper) J. Mater. Chem. A, 2019, 7, 312-321.

## Contributors
Jarvist Moore Frost and Suzanne K. Wallace. Codes started by Jarvist Moore Frost, September 19th 2014. Based on prior [Amphisbaena](https://github.com/jarvist/Amphisbaena) and [Starrynight](https://github.com/WMD-group/StarryNight) codes.
