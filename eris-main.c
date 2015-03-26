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

// Prototypes...

static double site_energy(int x, int y, int z, int species);
static void MC_move();

static int rand_int(int SPAN);

#include "eris-config.c" //Global variables & config file reader function  
#include "eris-lattice.c" //Lattice initialisation / zeroing / sphere picker fn; dot product
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

    initialise_lattice_random(); //populate wiht random dipoles


    fprintf(stderr,"Lattice initialised.\n");

    outputlattice_dumb_terminal(); //Party like it's 1980


    fprintf(stderr,"\n\tMC startup. 'Do I dare disturb the universe?'\n");

    fprintf(stderr,"'.' is %d MC moves attempted.\n",MCMinorSteps);

    fprintf(log,"# ACCEPT+REJECT, Efield, Eangle, E_dipole, E_strain, E_field, (E_dipole+E_strain+E_field)\n");

    //old code - now read in option, so I can parallise externally
    //    for (Efield.x=0.1; Efield.x<3.0; Efield.x+=0.5)
    //    for (T=0;T<1500;T+=100) //I know, I know... shouldn't hard code this.
    {
        beta=1/((float)T/300.0);

        for (i=0;i<TempSteps;i++)
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
            T=i*50;
            beta=1/((float)T/300.0);

            // Do some MC moves!

//            initialise_lattice_random();
            initialise_lattice_stripe();
            //#pragma omp parallel for //SEGFAULTS :) - non threadsafe code everywhere
            tic=time(NULL);
            for (j=0;j<MCMegaSteps;j++)
            {
                for (k=0;k<MCMinorSteps;k++) //let's hope the compiler inlines this to avoid stack abuse. Alternatively move core loop to MC_move fn?
                    MC_move();
            //    outputlattice_dumb_terminal();
            }
            toc=time(NULL);
 
            outputlattice_dumb_terminal(); //Party like it's 1980

            fprintf(stderr,"Efield: x %f y %f z %f | Dipole %f CageStrain %f K %f\n",Efield.x,Efield.y,Efield.z,Dipole,CageStrain,K);
            fflush(stdout); // flush the output buffer, so we can live-graph / it's saved if we interupt
            fprintf(stderr,"MC Moves: %f MHz\n",1e-6*(double)(MCMinorSteps*X*Y*Z)/(double)(toc-tic));
           
		lattice_potential_XYZ("potential.dat");
 
            // Manipulate the run conditions depending on simulation time
            //        if (i==100) { DIM=3;}  // ESCAPE FROM FLATLAND
            //        if (i==200) { Efield.z=1.0;}      // relax back to nothing
            //        if (i==300) {Efield.z=0.0; Efield.x=1.0;}

       }

    } 
    // OK; we're finished...

    fprintf(stderr,"\n");

    fprintf(stderr,"Monte Carlo moves - ACCEPT: %lu REJECT: %lu ratio: %f\n",ACCEPT,REJECT,(float)ACCEPT/(float)(REJECT+ACCEPT));
    fprintf(stderr," For us, there is only the trying. The rest is not our business. ~T.S.Eliot\n\n");

    return 0;
}

static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
    return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

static double site_energy(int x, int y, int z, int species)
{
    int dx,dy,dz=0;
    float d;
    double dE=0.0;

    int species2;

    // Sum over near neighbours for formalcharge-formalcharge interaction
    for (dx=-DipoleCutOff;dx<=DipoleCutOff;dx++)
        for (dy=-DipoleCutOff;dy<=DipoleCutOff;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-DipoleCutOff;dz<=DipoleCutOff;dz++) //NB: conditional zDipoleCutOff to allow for 2D version
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                d=sqrt((float) dx*dx + dy*dy + dz*dz); //that old chestnut; distance in Euler space
                if (d>(float)DipoleCutOff) continue; // Cutoff in d

                species2= lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

                dE+=E_int[species][species2]/d;
            }

    // Interaction of dipole with (unshielded) E-field
 /*   dE+= + dot(newdipole, & Efield)
        - dot(olddipole, & Efield);*/

    return(dE); 
}

static void MC_move()
{
    int x_a, y_a, z_a;
    int x_b, y_b, z_b;

    int dx, dy, dz;
    float d;
    float dE=0.0;

    int species_a, species_b;
    
    // Choose random dipole / lattice location

    x_a=rand_int(X);
    y_a=rand_int(Y);
    z_a=rand_int(Z);

    dx=(rand_int(2)*2)-1; // on interval [-1,1]
    dy=(rand_int(2)*2)-1;
    dz=(rand_int(2)*2)-1;

    x_b=(x_a+dx+X)%X;
    y_b=(y_a+dy+Y)%Y;
    z_b=(z_a+dz+Z)%Z;

    //if (lattice[x][y][z].x==0.0 && lattice[x][y][z].y==0.0 && lattice[x][y][z].z==0.0) return; //dipole zero length .'. not present

    // random new orientation. 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well

    species_a=lattice[x_a][y_a][z_a];
    species_b=lattice[x_b][y_b][z_b];

    //calc site energy
    // TODO: Check this! Self interaction? Species A vs. B? Want two
    // configuration states and diff in energy between them.
    dE+=site_energy(x_a,y_a,z_a, species_a);
    dE-=site_energy(x_a,y_a,z_a, species_b);

    dE+=site_energy(x_b,y_b,z_b, species_b);
    dE-=site_energy(x_b,y_b,z_b, species_a);

    // OK; now we have dE - proceed with Metropolis Accept/Reject Criteria
    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x_a][y_a][z_a]=species_b; //swap two atoms / species
        lattice[x_b][y_b][z_b]=species_a; 

        ACCEPT++;
    }
    else
        REJECT++;
}

