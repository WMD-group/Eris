/* Eris - a Monte Carlo code to simulate cation disorder in CZTS (Cu2ZnSnS4)
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

#include "mt19937ar-cok.c" //Code _included_ to allow more global optimisation
static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
        return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

#include "eris-config.c" //Global variables & config file reader function  
#include "eris-lattice.c" //Lattice initialisation / zeroing / sphere picker fn; dot product
#include "eris-kernel.c" // MC kernel, Lattice energies etc.
#include "eris-analysis.c" //Analysis functions, and output routines


// Functions for lattice initialisation and analysis of initial config

// Prototype functions
void initialise_lattice();

// wrapper function to put all logic in choosing lattice together in one place
void initialise_lattice()
{
    if (OrderedInitialLattice) // set by eris.cfg
    {
        if (SuzySupercell)
            initialise_lattice_CZTS_supercell();
        if (Method2)
            initialise_lattice_CZTS_method2();
        else
            initialise_lattice_CZTS();
    }
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

    sprintf(name,"Gulp_check_initial_T_%04d.in",T);
    if (EquilibrationChecksTestCutOff) generate_gulp_input(T, name);
    sprintf(name,"Eris_check_initial_all_r_T_%04d.dat",T);
    if (EquilibrationChecksTestCutOff) lattice_potential_r_test(name);

    sprintf(name,"T_%04d_gulp_input_initial.in",T);
    if (SaveGULP) generate_gulp_input(T, name);

    if (CalculateRadialOrderParameter) radial_distribution_function_allsites(0);
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
    if (ExtraData) generate_gulp_input(T, name);

    // Analysis and output routines
    if (DisplayDumbTerminal) outputlattice_dumb_terminal();
    if (EquilibrationChecks) { report_dE(); reset_dE(); }
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
        sprintf(filename,"GULP_inputs/T_%04d_gulp_input_final.in",T);
        generate_gulp_input(T, filename);
    }

    if (EquilibrationChecksTestCutOff)
    {
        char gulp_filename[100];
        sprintf(gulp_filename,"equilibration_check_Sn_potentials/Gulp_check_T_%04d.in",T);
        generate_gulp_input(T, gulp_filename);
        char eris_filename[100];
        sprintf(eris_filename,"equilibration_check_Sn_potentials/Eris_check_all_r_T_%04d.dat",T);
        lattice_potential_r_test(eris_filename);
    }

    // Output RDF of final configuration at each T
//    radial_distribution_function_allsites();
}



// Begin main loops of Monte Carlo simulation (calling above functions for analysis)

int main(int argc, char *argv[])
{
    int i,j,k, x,y; //for loop iterators
    int tic,toc,tac; // time counters; Goes the old grandfather clock

    fprintf(stderr,"Eris - Goddess of Kesterite Chaos.\n");

    fprintf(stderr,"Loading config...\n");
    load_config(); // see eris-config.c ; mainly sets global variables and Bools

    // Now override with command line options if supplied...
    // Currently assumes first argument is overriding TEMPERATURE; to be able
    // to be used with GNU parallel
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line temperature: Overiding settings and doing single shot calculation with: T = %d\n",T);
        TMAX=T; TMIN=T; // overide any loop settings from config; just this single temperature
    }

    // Allocate lattice; which is made up of 'dipole' structs
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

        {
            if (ReinitialiseLattice) 
            {
                initialise_lattice();
                fprintf(stderr,"Lattice reinitialised: \n");
            }
            else
                fprintf(stderr,"Lattice carried over from previous simulation: \n");

            analysis_initial();

            // Run the requested 'equilibration' steps as a burn in, before starting to
            // collect statistics 
                fprintf(stderr,"Equilibration Monte Carlo... (no data ouput): ");
                for (j=0;j<MCEqmSteps;j++)
                {
                    for (k=0;k<MCMinorSteps;k++)
                        MC_move();
                    fprintf(stderr,",");
                }
                fprintf(stderr,"\n");

            // START OF OUTER MC LOOP (MCMegaSteps)
              for (j=1;j<=MCMegaSteps;j++)
              {
                  // INNER MC LOOP (MCMinorSteps)
                  tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                  for (k=0;k<MCMinorSteps;k++) // Compiler should optimise this for loop out.
                      MC_move();
                  toc=clock();

                  analysis_midpoint(j); // analysis run after every MCMinorSteps block of moves 
                  fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output

                  tac=clock(); // timings for analysis/output

                  fprintf(stderr,"MC Moves: %f MHz\n",
                          1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC);
                  fprintf(stderr,"Time spent doing Monte Carlo as total fraction of time: %.2f %%\n",100.0*(double)(toc-tic)/(double)(tac-tic));
                  fprintf(stderr,"Monte Carlo moves - ATTEMPT: %llu ACCEPT: %llu REJECT: %llu Accept Ratio: %f\n",MCMinorSteps,ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
                  REJECT=0; ACCEPT=0;

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

