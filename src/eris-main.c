/* Eris: an on-lattice Monte Carlo code to simulate thermodynamic Cu-Zn disorder in kesterite-structured Cu2ZnSnS4 
 * 
 * Adapted from Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * Contributors: 
 * Suzanne Wallace, University of Bath
 *
 * File begun 16th January 2014
 */

#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h>

#include <sys/stat.h> // to enable to use of mkdir in functions in program
#include <sys/types.h>

// This code is _included_ (rather than compiled to objects) to allow more
// unrolling / function fusing.
#include "eris-random.c"   // PRNG
#include "eris-config.c"   // Global variables & config file reader function  
#include "eris-lattice.c"  // Lattice initialisation / zeroing / sphere picker fn; dot product
#include "eris-kernel.c"   // MC kernel, Lattice energies etc.
#include "eris-analysis.c" // Analysis functions, and output routines

// Functions for lattice initialisation and analysis of initial config

// Function prototypes
void initialise_lattice();
void analysis_initial();
void analysis_midpoint(int MCStep);
void analysis_final();
int main(int argc, char *argv[]);

// wrapper function to put all logic in choosing lattice together in one place
void initialise_lattice()
{
    if (OrderedInitialLattice) // set by eris.cfg
        initialise_lattice_CZTS_SKW();
    else
        initialise_lattice_CZTS_randomized();
}

// Analysis of the _initial_ lattice, before performing any MC moves 
void analysis_initial()
{
    char name[100];
    
    outputlattice_stoichometry(); // print histogram of stoichs for user; check to see what we have
    if (DisplayDumbTerminal) outputlattice_dumb_terminal(); // initial lattice

    // Producing GULP input file of lattice before performing MC moves for each T
    sprintf(name,"Sn_potentials_Temp_%04d_initial.dat",T); // for electrostatic potential file during equilibration check run
    if (EquilibrationChecks) equil_lattice_potential(name);

    sprintf(name,"Gulp_T_%04d_initial_lattice.in",T);
    if (SaveGULP) generate_gulp_input(T, name);

    if (CalculateRadialOrderParameter) radial_distribution_function_allsites_initial();
    if (SaveXYZ) outputlattice_xyz("czts_lattice_initial.xyz");
    if (CalculatePotential) lattice_potential_XYZ("potential_initial.dat");
//    lattice_energy(); // check energy sums
}


// Defining functions for analyses: during the simulation after each MC mega step (analysis_midpoint) and after each temperature step (analysis_final)

// OK, we've just done MCMinorSteps, now we can start analysis 
void analysis_midpoint(int MCStep)
{
    // Define log file names (stored as strings), with current temperature as part
    // of the filename
    char name[100];
    
    sprintf(name,"Sn_Pot_T_%04d_MCS_%05d.dat",T,MCStep); // for electrostatic potential file during equilibration check run
    if (EquilibrationChecks) equil_lattice_potential(name); 

    sprintf(name,"Gulp_T_%04d_MCS_%05d.in",T, MCStep);
    if (SaveGULP) generate_gulp_input(T, name);

    // Analysis and output routines
    if (DisplayDumbTerminal) outputlattice_dumb_terminal();
//    if (EquilibrationChecks) { report_dE(); reset_dE(); }  // Redundant now?
    if (CalculateRadialOrderParameter) radial_distribution_function_allsites(MCStep);

    sprintf(name,"potential_T_%04d_MCS_%05d.dat",T,MCStep); // for electrostatic potential file
    if (CalculatePotential) lattice_potential_XYZ(name);
}

// OK, we have now finished all of our MC steps for this T value
void analysis_final()
{
    if (SaveXYZ)
    {
        char name[100];
        sprintf(name,"czts_lattice_T_%04d.xyz",T);
        outputlattice_xyz(name);
    }

    if (SaveGULP)
    {
        char filename[100];
        sprintf(filename,"T_%04d_gulp_final_lattice.in",T);
        generate_gulp_input(T, filename);
    }

    if (SavePOSCAR)
    {
        char POSCARname[100];
        sprintf(POSCARname,"T_%04d_final_lattice.POSCAR.vasp",T);
        generate_POSCAR(POSCARname);
    }

    if (PotentialCubeFile)
    {
        char CubeFileName[100];
        sprintf(CubeFileName, "PotentialCubeFile_T_%04dK.cube", T);
        potential_3D_cube_file(CubeFileName);
    }

    if (BandTailingPotentials)
    {
        // Output Cu and Sn potentials of whole system at each T
        char CuPotFileName[100];
        char SnPotFileName[100];
        sprintf(CuPotFileName, "Cu_potentials_T_%04dK.dat", T);
        sprintf(SnPotFileName, "Sn_potentials_T_%04dK.dat", T);
        output_Cu_Sn_potentials(CuPotFileName, SnPotFileName);

        // Output potentials in 2D slices: Cu from Cu-Zn planes (odd slice numbers), Cu and Sn from Cu-Sn planes (even slice numbers)
        int z;
        for (z=0;z<Z;z++)
        {
            // Write files for potentials of Cu's in each Cu-Zn slice
            if (z%2 == 1) //for all odd z (i.e. only Cu-Zn layers from Method2 lattice initialisation)
            {
                char CuZnSlice[100]; 
                sprintf(CuZnSlice, "Cu_potentials_T_%04dK_slice_z=%02d.dat", T, z);
                CuZn_slice_potentials(CuZnSlice, z);
            }
            // Write files for potentials of Cu's and Sn's in each Cu-Sn slice
            if (z%2 == 0) //for even odd z (i.e. only Cu-Sn layers from Method2 lattice initialisation)
            {
                char CuSnSlice_Cu[100]; 
                sprintf(CuSnSlice_Cu, "Cu_potentials_T_%04dK_slice_z=%02d.dat", T, z);
                char CuSnSlice_Sn[100]; 
                sprintf(CuSnSlice_Sn, "Sn_potentials_T_%04dK_slice_z=%02d.dat", T, z);
                CuSn_slice_potentials(CuSnSlice_Cu, CuSnSlice_Sn, z);
            }
        }
    }



    // Output RDF of final configuration at each T
//    radial_distribution_function_allsites();
}


// Main function:-
// * Read in config from eris.cfg
// * Override with command line options if specified (argc,argv etc.)
// * Run main MC loop, running the analysis functions above
// * Exit once finished
int main(int argc, char *argv[])
{
    int i,j,k, x,y; //for loop iterators
    int tic,toc,tac; // time counters; Goes the old grandfather clock

    fprintf(stderr,"Eris - Goddess of Kesterite Chaos.\n");

    fprintf(stderr,"Loading config...\n");
    load_config(); // see eris-config.c ; sets global variables and Bools

    // Now override with command line options if supplied...
    // Assumes *first* command line argument is overriding TEMPERATURE; 
    // mainly for (non-interactive) use with GNU parallel
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line temperature: Overiding settings and doing single shot calculation with: T = %d\n",T);
        TMAX=T; TMIN=T; // overide any loop settings from config; just this single temperature
    }

    // Allocate lattice
    fprintf(stderr,"Memory allocation for lattice with X=%d Y=%d Z=%d\n",X,Y,Z);
    lattice = (int ***)malloc(sizeof(int **)*X);
    for (x=0;x<X;x++) 
    {
        lattice[x]=(int **)malloc(sizeof(int *)*Y);
        for (y=0;y<Y;y++)
            lattice[x][y]=(int *)malloc(sizeof(int)*Z);
    } 
    fprintf(stderr,"Lattice allocated");

    // If we're going to do some actual science, we better have a logfile...
    FILE *log;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    // TODO: Make this a setting in eris.cfg.
    //Fire up the Mersenne twister!
    int SEED=0xDEADBEEF;
    //int SEED=time(NULL);
    init_genrand(SEED); // Initialise Pseudo-Random-Number-Generator
    fprintf(stderr,"Mersenne Twister initialised with seed: %X. \n",SEED);
    fprintf(log,"# Mersenne Twiser initialised with seed: %X",SEED);
    fprintf(log,"# Eris run started at time = %ld",time(NULL));

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");
    fprintf(stderr,"'.' is %llu MC moves attempted.\n",MCMinorSteps);

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    initialise_lattice(); // contains logic about choosing what lattice to use   

    analysis_initial();

    // Core simulation loop
    for (T=TMAX;T>=TMIN;T-=TSTEP) // read in from eris.cfg; or overidden by single shot T from command line 
    {
        beta=1/((float)T/300.0); // Thermodynamic Beta; in units of kbT @ 300K
        printf("Temperature now T: %d K \t Beta: %f (kbT@300K)\n",T,beta);

            if (ReinitialiseLattice) 
            {
                initialise_lattice();
                fprintf(stderr,"Lattice reinitialised: \n");
            }
            else
                fprintf(stderr,"Lattice carried over from previous simulation: \n");

            analysis_initial();

            // Run the requested 'equilibration' steps as a burn in, before starting to collect statistics
            if (!(EquilibrationChecks)) // Do not perform if doing an equilibration test
            { 
                fprintf(stderr,"Equilibration Monte Carlo... (no data ouput): ");
                for (j=0;j<MCEqmSteps;j++)
                {
                    for (k=0;k<MCMinorSteps;k++)
                        MC_move_stencil();
                    fprintf(stderr,",");
                }
                fprintf(stderr,"\n");
            }
            // Perform electrostatics check instead of standard MC simulation if this flag is set in eris.cfg
            if (ElectrostaticsCheck)
            {
                for (j=1;j<=MCMegaSteps;j++)
                {
                    // INNER MC LOOP (MCMinorSteps)
                    tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                    //for (k=0;k<MCMinorSteps;k++) // Compiler should optimise this for loop out. //Remove for dE check - too many moves outputted!
                    MC_move_dE_check_stencil();
                    toc=clock();
                    fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output
                    tac=clock(); // timings for analysis/output
                    fprintf(stderr,"Performing electrostatics check...\n");
                    fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
                }
            }
            else
            {
            // START OF OUTER MC LOOP (MCMegaSteps)
                for (j=1;j<=MCMegaSteps;j++)
                {
                    // INNER MC LOOP (MCMinorSteps)
                    tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                    for (k=0;k<MCMinorSteps;k++) // Compiler should optimise this for loop out.
                        MC_move_stencil();
                    toc=clock();
  
                    analysis_midpoint(j); // analysis run after every MCMinorSteps block of moves 
                    fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output
  
                    tac=clock(); // timings for analysis/output
                    fprintf(stderr,"MC Moves: %f MHz\n",
                    1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC);
                    fprintf(stderr,"Time spent doing Monte Carlo as total fraction of time: %.2f %%\n",100.0*(double)(toc-tic)/(double)(tac-tic));
                    fprintf(stderr,"Monte Carlo moves - ATTEMPT: %llu ACCEPT: %llu REJECT: %llu Accept Ratio: %f\n",MCMinorSteps,ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
                    REJECT=0; ACCEPT=0; // reset counters
  
                    fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
                }
                analysis_final(); //analysis run at end of sim
            }
        }

    // OK; we're finished...
    fprintf(stderr,"\n");
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

