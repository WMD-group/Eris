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

int main(int argc, char *argv[])
{
    int i,j,k, x,y; //for loop iterators
    int tic,toc,tac; // time counters; Goes the old grandfather clock
    
    double P=0.0;

//---------------------------------------------------------------------------------------------------------------------------------- 
// Initial terminal messages to user, reading in config files and initialising CZTS lattice
//---------------------------------------------------------------------------------------------------------------------------------- 

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

    // If we're going to do some actual science, we better have a logfile...
    FILE *log;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    //Fire up the Mersenne twister!
    init_genrand(0xDEADBEEF); //314159265);  // reproducible data :)
    //init_genrand(time(NULL)); // seeded with current time
    fprintf(stderr,"Mersenne Twister initialised. ");

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");
    fprintf(stderr,"'.' is %llu MC moves attempted.\n",MCMinorSteps);

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    if (OrderedInitialLattice) // set by eris.cfg
    {
        if (SuzySupercell)
            initialise_lattice_CZTS_supercell();
        else
            initialise_lattice_CZTS();
    }
    else
        initialise_lattice_CZTS_randomized();
    
//    outputlattice_stoichometry(); // print histogram of stoichs for user; check to see what we have

    if (DisplayDumbTerminal) outputlattice_dumb_terminal(); // initial lattice

    // Analysis files of the initial lattice, before performing any MC moves 
    if (CalculateRadialOrderParameter) radial_distribution_function_allsites_initial();
    if (SaveXYZ) outputlattice_xyz("czts_lattice_initial.xyz");
    if (CalculatePotential) lattice_potential_XYZ("potential_initial.dat");
	
//    lattice_energy(); // check energy sums


//---------------------------------------------------------------------------------------------------------------------------------- 
// START OF OUTERMOST TEMPERATURE LOOP
//---------------------------------------------------------------------------------------------------------------------------------- 

    // Core simulation loop
    for (T=TMAX;T>=TMIN;T-=TSTEP) // read in from eris.cfg 
    {
        beta=1/((float)T/300.0); // Thermodynamic Beta; in units of kbT @ 300K
        printf("Temperature now T: %d K \t Beta: %f (kbT@300K)\n",T,beta);

        {

// Re-initialising lattice for new temperature (if requested by user)
//---------------------------------------------------------------------------------------------------------------------------------- 
            if (ReinitialiseLattice) // Are we intending to reset the lattice?
            {
                if (OrderedInitialLattice)
                {
                    if (SuzySupercell)
                        initialise_lattice_CZTS_supercell();
                    else
                        initialise_lattice_CZTS();
                }
                else
                    initialise_lattice_CZTS_randomized();

                fprintf(stderr,"Lattice reinitialised: \n");
            }
            else
                fprintf(stderr,"Lattice carried over from previous simulation: \n");
           
            if (DisplayDumbTerminal) outputlattice_dumb_terminal();


// Defining file names based on temperature
//---------------------------------------------------------------------------------------------------------------------------------- 

            char electrostaticpotential_filename[100];
            sprintf(electrostaticpotential_filename,"potential_T_%04d.dat",T); // for electrostatic potential file
            char electrostaticpotential_equil_filename[100];
            sprintf(electrostaticpotential_equil_filename,"equilibration_check_potential+variance/equil_potential_T_%04d.dat",T); // for electrostatic potential file during equilibration check run
            char variance_equil_filename[100];
            sprintf(variance_equil_filename,"equilibration_check_potential+variance/equil_variance_T_%04d.dat",T); // for variance of potential file during equilibriation run as a function of MC step (or j in MCMegaSteps loop)

            // Producing GULP input file of lattice before performing MC moves for each T
            if (EquilibrationChecks) 
            {
                char gulp_filename_initial[100];
                sprintf(gulp_filename_initial,"equilibration_check_GULP_inputs/gulp_input_Temp_%04d_initial.in",T);
                lattice_energy_full(gulp_filename_initial);
            }

            // Producing gulp input file of initial lattice
            if (SaveGULP)
            {
                char filename[100];
                sprintf(filename,"GULP_inputs/czts_lattice_initial_T_%04d.in",T);
                generate_gulp_input(filename);
            }

// Jarv's equilibration run to use after suitable equilibration time has been established, set using MCEqmSteps in eris.cfg
//---------------------------------------------------------------------------------------------------------------------------------- 

            // equilibration before data collection
            if (MCEqmSteps>0)
            {
                fprintf(stderr,"Equilibration Monte Carlo... (no data ouput): ");
                for (j=0;j<MCEqmSteps;j++)
                {
                    for (k=0;k<MCMinorSteps;k++)
                        MC_move();
                    fprintf(stderr,",");
                }
                fprintf(stderr,"\n");
            }

//---------------------------------------------------------------------------------------------------------------------------------- 
// START OF OUTER MC LOOP (MCMegaSteps)
//---------------------------------------------------------------------------------------------------------------------------------- 
           
            //#pragma omp parallel for //SEGFAULTS :) - non threadsafe code everywhere
            for (j=0;j<MCMegaSteps;j++)
            {
//---------------------------------------------------------------------------------------------------------------------------------- 
// INNER MC LOOP (MCMinorSteps)
//---------------------------------------------------------------------------------------------------------------------------------- 
                tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                for (k=0;k<MCMinorSteps;k++) // Compiler should optimise this for loop out.
                    MC_move();
//---------------------------------------------------------------------------------------------------------------------------------- 
                toc=clock();


// Equilibration checks, performed for each MCMegaStep, after each complete loop of MCMinorSteps
//---------------------------------------------------------------------------------------------------------------------------------- 
                if (EquilibrationChecks) 
                {
                   T_separated_lattice_potential(electrostaticpotential_equil_filename, variance_equil_filename, j);
                   // Generating gulp input files for intermittent configurations during equilibriation for post-processing to calculate full lattice energy with gulp
                   char gulp_filename[100];
                   sprintf(gulp_filename,"equilibration_check_GULP_inputs/gulp_input_Temp_%04d_MCS_%04d.in",T,j);
                   lattice_energy_full(gulp_filename);
                 }
//---------------------------------------------------------------------------------------------------------------------------------- 

                
// Analysis and output routines
//---------------------------------------------------------------------------------------------------------------------------------- 
                if (DisplayDumbTerminal) outputlattice_dumb_terminal();
                if (CalculateRadialOrderParameter) radial_distribution_function_allsites();
		        if (CalculatePotential) lattice_potential_XYZ(electrostaticpotential_filename);

                fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output

                tac=clock(); // timings for analysis/output

                fprintf(stderr,"MC Moves: %f MHz\n",
                    1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC);
                fprintf(stderr,"Time spent MC vs. analysis: %f %\n",100.0*(double)(toc-tic)/(double)(tac-tic));

                fprintf(stderr,"Monte Carlo moves - ATTEMPT: %llu ACCEPT: %llu REJECT: %llu Accept Ratio: %f\n",MCMinorSteps,ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
                REJECT=0; ACCEPT=0;

                fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            }

//---------------------------------------------------------------------------------------------------------------------------------- 
            // OK, we have now finished all of our MC steps for this T value
            if (SaveXYZ)
            {
                char name[100];
                sprintf(name,"czts_lattice_T_%04d.xyz",T);
                outputlattice_xyz(name);
            }
       
        
            if (SaveGULP)
            {
                char filename[100];
                sprintf(filename,"GULP_inputs/czts_lattice_final_T_%04d.in",T);
                generate_gulp_input(filename);
            }
        
            // Output RDF of final configuration at each T
            radial_distribution_function_allsites();

        }

    }

    // OK; we're finished...
    fprintf(stderr,"\n");
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

