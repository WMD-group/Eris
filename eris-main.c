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

    initialise_lattice_CZTS();
 
    //old code - now read in option, so I can parallise externally
    //    for (Efield.x=0.1; Efield.x<3.0; Efield.x+=0.5)
        for (T=10000;T>0;T=T*0.99) //I know, I know... shouldn't hard code this.
    {
        beta=1/((float)T/300.0);
        printf("T: %d beta: %f\n",T,beta);

//        for (i=0;i<TempSteps;i++)
        {
            // Alright, this is the plan
            // First we take our variable
            // Then we bit reverse it as binary
            // Pretty confusing, but means it will fill in the temperature
            // range with maximum coverage, rather than linear ramping
          /*  unsigned char r=i;
            r=(r&0xF0)>>4 | (r&0x0F)<<4;
            r=(r&0xCC)>>2 | (r&0x33)<<2;
            r=(r&0xAA)>>1 | (r&0x55)<<1;
*/
//            T=i*50;
//            beta=1/((float)T/300.0);

            // Do some MC moves!

            char electrostaticpotential_filename[100];
            sprintf(electrostaticpotential_filename,"potential_T_%04d.dat",T); // for electrostatic potential file

//            initialise_lattice_random();
//            initialise_lattice_stripe();
//            initialise_lattice_CZTS();
            // test RDF routine...
            radial_distribution_function();
            fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output

            fprintf(stderr,"Lattice initialised.\n");
            outputlattice_xyz("czts_lattice_initial.xyz");
	        outputlattice_dumb_terminal();
//break;
            //#pragma omp parallel for //SEGFAULTS :) - non threadsafe code everywhere
            for (j=0;j<MCMegaSteps;j++)
            {
                tic=clock(); // measured in CLOCKS_PER_SECs of a second. 
                for (k=0;k<MCMinorSteps;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
                    MC_move();
                toc=clock();

// Analysis and output routines
                outputlattice_dumb_terminal();
//                radial_distribution_function();
                fflush(stdout); // flush buffer, so data is pushed out & you can 'ctrl-c' the program, retaining output
                fprintf(stderr,"MC Moves: %f MHz\n",
                    1e-6*(double)(MCMinorSteps)/(double)(toc-tic)*(double)CLOCKS_PER_SEC); 
                fprintf(stderr,"Monte Carlo moves - ATTEMPT: %llu ACCEPT: %llu REJECT: %llu Accept Ratio: %f\n",MCMinorSteps,ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
                REJECT=0; ACCEPT=0;

                fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            }
 
            fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
		    lattice_potential_XYZ(electrostaticpotential_filename);

            char name[100];
            sprintf(name,"czts_lattice_T_%04d.xyz",T);
            outputlattice_xyz(name);
            // Manipulate the run conditions depending on simulation time
            //        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
            //        if (i==200) { Efield.z=1.0;}      // relax back to nothing
            //        if (i==300) {Efield.z=0.0; Efield.x=1.0;}
    
       }

    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    fprintf(stderr,"Monte Carlo moves - ACCEPT: %llu REJECT: %llu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

