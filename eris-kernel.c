/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

static double site_energy(int x, int y, int z, int species, int CutOff);
static void MC_move();

static void equlibriation_statistics(float dE);
double sum_dE=0.0;

static double site_energy(int x, int y, int z, int species_a, int CutOff)
{
    int dx,dy,dz=0;
    float d;
    double dE=0.0;

    int species_b;

    // Sum over near neighbours for formalcharge-formalcharge interaction
#pragma omp parallel for reduction(+:dE) 
// OPENMP PARALLISATION
    for (dx=-CutOff;dx<=CutOff;dx++)
        for (dy=-CutOff;dy<=CutOff;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-CutOff;dz<=CutOff;dz++) //NB: conditional CutOff to allow for 2D version
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                d=sqrt((float) dx*dx + dy*dy + dz*dz); //that old chestnut; distance in Euler space
//                if (d>(float)DipoleCutOff) continue; // Cutoff in d
//                -->
//                EXPANSIONS IN SPHERES; probably not convergent

                species_b= lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

                // E_int runs from 1..4
                if (species_a>0 && species_b>0) // if gaps in lattice, no interaction energy contribution
                    dE+=E_int[species_a-1][species_b-1]/d;
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
    
    // Choose random _occupied_ lattice location
    do 
    {
        x_a=rand_int(X);
        y_a=rand_int(Y);
        z_a=rand_int(Z);
    }
    while (lattice[x_a][y_a][z_a]==0); // keep selecting new random numbers until find occupied site...

    // Choice direction + size of moves...
    int RADIAL_CUTOFF=3;
    do
    {
// Local cube; Anyhere <> limit
//        dx=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
//        dy=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
//        dz=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;

// Nearest neighbour over
//    dx=(rand_int(2)*2)-1; // on interval [-1,1]
//    dy=(rand_int(2)*2)-1;
//    dz=(rand_int(2)*2)-1;

// Global move; anywhere <> anywhere in lattice
    dx=rand_int(X);
    dy=rand_int(Y);
    dz=rand_int(Z);
  
    }
    while( (dx+dy+dz)%2!=0 || dx==dy==dz==0); //check to see whether site at this offset in gappy FCC lattice.
    // and we're not trying to swap with ourselves...

//    fprintf(stderr,"MC_move: dx dy dz %d %d %d\n",dx,dy,dz); 
        // for debug - check weighting of moves

//    if (dx==dy==dz==0) return; //skip consideration if null move

    // 2nd site to look at...
    x_b=(x_a+dx+X)%X;
    y_b=(y_a+dy+Y)%Y;
    z_b=(z_a+dz+Z)%Z;

    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well

    species_a=lattice[x_a][y_a][z_a];
    species_b=lattice[x_b][y_b][z_b];

    // Immobilise Tin!
    // Nb: Not very computationally efficient, better would be to never select
    // Tin in the first place...
//    if (species_a==4 || species_b==4) // if either move selects Tin...
//        return;

    if (species_a==0 || species_b==0) // if interstial / empty site...
        return; // don't do a move. Highly computational inefficient, FIXME

    if (species_a==species_b) // move achieves nothing... don't count both
       // calculating NULL moves, or counting them towards ACCEPT/REJECT criter
            return;

    //calc site energy
    // TODO: Check this! Self interaction? Species A vs. B? Want two
    // configuration states and diff in energy between them.
    dE+=site_energy(x_a,y_a,z_a, species_a, ElectrostaticCutOff);
    dE-=site_energy(x_a,y_a,z_a, species_b, ElectrostaticCutOff);

    dE+=site_energy(x_b,y_b,z_b, species_b, ElectrostaticCutOff);
    dE-=site_energy(x_b,y_b,z_b, species_a, ElectrostaticCutOff);

    // Report on planned move + dE -- for debugging only (makes a ridiculous
    // number of prints...)
/*    fprintf(stderr,"MC Move: %d on %d %d %d, %d on %d %d %d --> dE: %f\n",
            lattice[x_a][y_a][z_a], x_a, y_a, z_a,
            lattice[x_b][y_b][z_b], x_b, y_b, z_b,
            dE);
*/
    // OK; now we have dE - proceed with Metropolis Accept/Reject Criteria
    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        lattice[x_a][y_a][z_a]=species_b; //swap two atoms / species
        lattice[x_b][y_b][z_b]=species_a; 

        ACCEPT++;

//        equlibriation_statistics(dE); # Pausing dev on this while making new
//        array structures
    }
    else
        REJECT++;
}

#define BINS 100
int endo_bins[BINS], exo_bins[BINS];

static void equlibriation_statistics(float dE)
{
     sum_dE+=dE;

     if (ACCEPT%10000==0) // every %XXX accepted moves...
     {
         printf("Change in energy over 1E4 moves: %f\n",sum_dE);
         sum_dE=0.0;
     }


     if (dE>0.0)
     {
//     printf("beta: %f dE: %f logf(fabsf(dE)): %f logf(fabsf(dE*beta)): %f\n",beta,dE,logf(fabsf(dE)),pow(fabsf(dE),0.5));
        endo_bins[ (int) (pow(fabsf(dE),0.5)/0.05) ] ++; 
     }
     else
         exo_bins[ (int) (pow(fabsf(dE),0.5)/0.05) ] ++;

     if (ACCEPT%100000==0)
     {
        printf("Histogram of exo/endo over 1E5 moves: \n");
         for (int i=0; i<BINS; i++)
         {
             printf("%d",endo_bins[i]);
             printf("/%d ",exo_bins[i]);
             endo_bins[i]=0;
             exo_bins[i]=0;
         }
     }
}

