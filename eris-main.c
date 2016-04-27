/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <math.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <libconfig.h>

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

    double P=0.0;

    int tic,toc;

    fprintf(stderr,"Eris - Goddess of Kesterite Chaos.\n");

    fprintf(stderr,"Loading config...\n");
    load_config();

    // Now override with command line options if supplied...
    if (argc>1)
    {
        sscanf(argv[1],"%d",&T);
        fprintf(stderr,"Command line temperature: T = %d\n",T);
    }

    // If we're going to do some actual science, we better have a logfile...
    FILE *log;
    log=fopen(LOGFILE,"w");
    fprintf(stderr,"Log file '%s' opened. ",LOGFILE);

    //Fire up the twister!
    init_genrand(0xDEADBEEF); //314159265);  // reproducible data :)
    //init_genrand(time(NULL)); // seeded with current time
    fprintf(stderr,"Twister initialised. ");

    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");
    fprintf(stderr,"'.' is %llu MC moves attempted.\n",MCMinorSteps);

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    if (OrderedInitialLattice)
        initialise_lattice_CZTS();
    else
        initialise_lattice_CZTS_randomized();
    
    outputlattice_stoichometry(); // print histogram of stoichs for user.

    if (DisplayDumbTerminal) outputlattice_dumb_terminal(); // initial lattice

    if (CalculateRadialOrderParameter) radial_distribution_function("RDF_initial.dat");

    lattice_energy(); // check energy sums
    //exit(-1);

    //old code - now read in option, so I can parallise externally
    //    for (Efield.x=0.1; Efield.x<3.0; Efield.x+=0.5)
    for (T=TMAX;T>=TMIN;T-=TSTEP) // read in from eris.cfg 
    {
        beta=1/((float)T/300.0);
        printf("T: %d beta: %f\n",T,beta);

        {
            // Do some MC moves!

            char electrostaticpotential_filename[100];
            sprintf(electrostaticpotential_filename,"potential_T_%04d.dat",T); // for electrostatic potential file
            // Setting filename for RDF output
            char RDF_filename[100];
            sprintf(RDF_filename,"RDF_T_%04d.dat",T); // for electrostatic potential file
            
            
            if (ReinitialiseLattice) // Are we intending to reset the lattice?
            {
                if (OrderedInitialLattice)
                    initialise_lattice_CZTS();
                else
                    initialise_lattice_CZTS_randomized();

                fprintf(stderr,"Lattice reinitialised at T = %d K\n",T);
            }
            else
                fprintf(stderr,"Lattice carried over from previous simulation. Now at T = %d K\n",T);

            if (CalculateRadialOrderParameter) radial_distribution_function(RDF_filename);

            //            fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output

            // Outputting xyz file for initial lattice and potential file before performing any Monte Carlo moves
            if (SaveXYZ) outputlattice_xyz("czts_lattice_initial.xyz");
            if (CalculatePotential) lattice_potential_XYZ("potential_initial.dat");
	        
            
//            if (DisplayDumbTerminal) outputlattice_dumb_terminal();
//break;
            //#pragma omp parallel for //SEGFAULTS :) - non threadsafe code everywhere
            for (j=0;j<MCMegaSteps;j++)
            {
                tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                for (k=0;k<MCMinorSteps;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
                    MC_move();
                toc=clock();

// Analysis and output routines
                if (DisplayDumbTerminal) outputlattice_dumb_terminal();
                if (CalculateRadialOrderParameter) radial_distribution_function(RDF_filename);

                fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output
                fprintf(stderr,"MC Moves: %f MHz\n",
                    1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC); 
                fprintf(stderr,"Monte Carlo moves - ATTEMPT: %llu ACCEPT: %llu REJECT: %llu Accept Ratio: %f\n",MCMinorSteps,ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
                REJECT=0; ACCEPT=0;

                fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            }
 
            fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
		    if (CalculatePotential) lattice_potential_XYZ(electrostaticpotential_filename);

            if (SaveXYZ)
            {
                char name[100];
                sprintf(name,"czts_lattice_T_%04d.xyz",T);
                outputlattice_xyz(name);
            }
       }

    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    fprintf(stderr,"Monte Carlo moves - ACCEPT: %llu REJECT: %llu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

