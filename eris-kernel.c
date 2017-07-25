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
static void MC_move_dE_check();
static void log_dE(float dE);
static int min(int a, int b, int c);



#define EVJEN 1 // could convert these into ints later if wanting to make them dynamic
#define SPHERICAL 0 
// Nb: in C, 0=FALSE, 1=TRUE



static int min(int a, int b, int c)
{
  int m = a;
  if (m > b) m = b;
  if (m > c) m = c;
  return m;
}





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

                if (SPHERICAL)
                    if (d>(float)CutOff) continue; // Cutoff in d
//                -->
//                EXPANSIONS IN SPHERES; probably not convergent

                species_b= lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];

                // E_int runs from 1..4
                if (species_a==0 || species_b==0) // if gaps in lattice, no interaction energy contribution
                    continue;

                // "... the potentials of of the ions forming the surface of
                // the cube, however, are given the weights 1/2, 1/4 or 1/8
                // according as they are situated on a face, an edge, or
                // a corner of the cube.
                // Evjen - Physical Review Vol 39, 1932
                double evjen_E=E_int[species_a-1][species_b-1]/d;
                double evjen_weight=1.0;
                if (EVJEN)
                {
                    if (abs(dx)==CutOff) evjen_weight*=0.5; 
                    if (abs(dy)==CutOff) evjen_weight*=0.5;
                    if (abs(dz)==CutOff) evjen_weight*=0.5;// corner
//                  fprintf(stderr,"X:%d Y:%d Z:%d CutOff:%d evjenweight: %f\n",dx,dy,dz,CutOff,evjen_weight);
                }
                dE+=evjen_E * evjen_weight;
            }

    // Interaction of dipole with (unshielded) E-field
 /*   dE+= + dot(newdipole, & Efield)
        - dot(olddipole, & Efield);*/

    return(dE); 
}

static double site_energy_stencil(int x, int y, int z, int species_a, int CutOff, int sx, int sy, int sz)
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
//                if (dx==0 && dy==0 && dz==0)
//                    continue; //no infinities / self interactions please!

//                d=sqrt((float) dx*dx + dy*dy + dz*dz); //that old chestnut; distance in Euler space

                d=sqrt( (float) (sx+dx-x)*(sx+dx-x) + (sy+dy-y)*(sy+dy-y) + (sz+dz-z)*(sz+dz-z) );
                
                if (d<0.5) continue;

//                if (SPHERICAL)
//                    if (d>(float)CutOff) continue; // Cutoff in d
//                -->
//                EXPANSIONS IN SPHERES; probably not convergent

                species_b= lattice[(X+sx+dx)%X][(Y+sy+dy)%Y][(Z+sz+dz)%Z];

                // E_int runs from 1..4
                if (species_a==0 || species_b==0) // if gaps in lattice, no interaction energy contribution
                    continue;

                // "... the potentials of of the ions forming the surface of
                // the cube, however, are given the weights 1/2, 1/4 or 1/8
                // according as they are situated on a face, an edge, or
                // a corner of the cube.
                // Evjen - Physical Review Vol 39, 1932
                double evjen_E=E_int[species_a-1][species_b-1]/d;
                double evjen_weight=1.0;
                if (EVJEN)
                {
                    if (abs(dx)==CutOff) evjen_weight*=0.5; 
                    if (abs(dy)==CutOff) evjen_weight*=0.5;
                    if (abs(dz)==CutOff) evjen_weight*=0.5;// corner
//                  fprintf(stderr,"X:%d Y:%d Z:%d CutOff:%d evjenweight: %f\n",dx,dy,dz,CutOff,evjen_weight);
                }
                dE+=evjen_E * evjen_weight;
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
    int RADIAL_CUTOFF=2;
    do
    {
// Local cube; Anyhere <> limit
// Nearest Neighbour limit is RADIAL_CUTOFF=1
        dx=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        dy=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        //dz=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        dz=0; // freeze motion in Z; i.e. between Cu/Zn and Cu/Sn layers
    }
    while( (dx==0 && dy==0 && dz==0) || (dx+dy+dz)%2!=0 ); // Check this works as intended!
        // check to see whether site at this offset in gappy FCC lattice.
    // and we're not trying to swap with ourselves...

    //fprintf(stderr,"MC_move: dx dy dz %d %d %d\n",dx,dy,dz);
        // for debug - check weighting of moves

    // 2nd site to look at...
    x_b=(x_a+dx+X)%X;
    y_b=(y_a+dy+Y)%Y;
    z_b=(z_a+dz+Z)%Z;

    // Global move; anywhere <> anywhere in lattice
    //x_b=rand_int(X);
    //y_b=rand_int(Y);
    //z_b=rand_int(Z);

    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well

    species_a=lattice[x_a][y_a][z_a];
    species_b=lattice[x_b][y_b][z_b];

    // Immobilise Tin!
    // Nb: Not very computationally efficient, better would be to never select
    // Tin in the first place...
    if (freezeSn)
        if (species_a==Sn || species_b==Sn) // if either move selects Tin...
            return;

    if (species_a==0 || species_b==0) // if interstial / empty site...
        return; // don't do a move. Highly computational inefficient, FIXME

    if (species_a==species_b) // move achieves nothing... don't count both
       // calculating NULL moves, or counting them towards ACCEPT/REJECT criter
            return;

    //calc site energy
    // TODO: Check this! Self interaction? Species A vs. B? Want two
    // configuration states and diff in energy between them.
//    dE+=site_energy(x_a,y_a,z_a, species_a, ElectrostaticCutOff);   
//     dE-=site_energy(x_a,y_a,z_a, species_b, ElectrostaticCutOff);

//    dE+=site_energy(x_b,y_b,z_b, species_b, ElectrostaticCutOff);
//    dE-=site_energy(x_b,y_b,z_b, species_a, ElectrostaticCutOff);

    dE+=site_energy_stencil(x_a,y_a,z_a, species_a, ElectrostaticCutOff, x_a, y_a, z_a);
    dE-=site_energy_stencil(x_a,y_a,z_a, species_b, ElectrostaticCutOff, x_a, y_a, z_a);

    dE+=site_energy_stencil(x_b,y_b,z_b, species_b, ElectrostaticCutOff, x_a, y_a, z_a);
    dE-=site_energy_stencil(x_b,y_b,z_b, species_a, ElectrostaticCutOff, x_a, y_a, z_a);


    // Report on planned move + dE -- for debugging only (makes a ridiculous
    // number of prints...)
    if (DEBUG)
        fprintf(stderr,"MC Move: %d on %d %d %d, %d on %d %d %d --> dE: %f\n",
                lattice[x_a][y_a][z_a], x_a, y_a, z_a,
                lattice[x_b][y_b][z_b], x_b, y_b, z_b,
                dE);

    // OK; now we have dE - proceed with Metropolis Accept/Reject Criteria
    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        //ACCEPT move
        lattice[x_a][y_a][z_a]=species_b; //swap two atoms / species
        lattice[x_b][y_b][z_b]=species_a;

        ACCEPT++;

        if (EquilibrationChecks) log_dE(dE);
    }
    else
        //REJECT move; doesn't need to do anything, just update REJECT for
        //stats
        REJECT++;
}






static void MC_move_dE_check()
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
    int RADIAL_CUTOFF=2;
    do
    {
// Local cube; Anyhere <> limit
// Nearest Neighbour limit is RADIAL_CUTOFF=1
        dx=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        dy=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        //dz=(rand_int(1+RADIAL_CUTOFF*2))-RADIAL_CUTOFF;
        dz=0; // freeze motion in Z; i.e. between Cu/Zn and Cu/Sn layers
    }
    while( (dx==0 && dy==0 && dz==0) || (dx+dy+dz)%2!=0 ); // Check this works as intended! 
        // check to see whether site at this offset in gappy FCC lattice.
    // and we're not trying to swap with ourselves...

    //fprintf(stderr,"MC_move: dx dy dz %d %d %d\n",dx,dy,dz); 
        // for debug - check weighting of moves


    // 2nd site to look at...
    x_b=(x_a+dx+X)%X;
    y_b=(y_a+dy+Y)%Y;
    z_b=(z_a+dz+Z)%Z;

    // Global move; anywhere <> anywhere in lattice
    //x_b=rand_int(X);
    //y_b=rand_int(Y);
    //z_b=rand_int(Z);
 
    // Nb: this is the definition of a MC move - might want to consider
    // alternative / global / less disruptive moves as well

    species_a=lattice[x_a][y_a][z_a];
    species_b=lattice[x_b][y_b][z_b];

    // Immobilise Tin!
    // Nb: Not very computationally efficient, better would be to never select
    // Tin in the first place...
    if (freezeSn)
        if (species_a==Sn || species_b==Sn) // if either move selects Tin...
            return;

    if (species_a==0 || species_b==0) // if interstial / empty site...
        return; // don't do a move. Highly computational inefficient, FIXME


    if (species_a==species_b) // move achieves nothing... don't count both
       // calculating NULL moves, or counting them towards ACCEPT/REJECT criter
          return;



// original site_energy_stencil to calculate dE ------------------------------------------------

//    dE+=site_energy_stencil(x_a,y_a,z_a, species_a, ElectrostaticCutOff, x_a, y_a, z_a);
//    dE-=site_energy_stencil(x_a,y_a,z_a, species_b, ElectrostaticCutOff, x_a, y_a, z_a);

//    dE+=site_energy_stencil(x_b,y_b,z_b, species_b, ElectrostaticCutOff, x_a, y_a, z_a);
//    dE-=site_energy_stencil(x_b,y_b,z_b, species_a, ElectrostaticCutOff, x_a, y_a, z_a);
// ----------------------------------------------------------------------------------------------

// Suzy add in: dE convergence check wrt ElectrostaticCutoff - add in loop increasing radius and printing out dE for run set to 1 MC step?

    FILE *fo;
    char dE_check_file[100];
    sprintf(dE_check_file,"dE_conv_check_T_%04d.dat",T);
    fo=fopen(dE_check_file,"a");
    fprintf(fo,"#r_cutoff dE \n");
    int CutOff_test=0;

    // Defining max for cut off radius for lattice sums based on lattice dimensions (i.e. min dimension/2) 
    int CutOffMax = floor(min(X,Y,Z)/2);

    for (CutOff_test=1; CutOff_test<CutOffMax+1; CutOff_test++)
    {
      dE+=site_energy(x_a,y_a,z_a, species_a, CutOff_test);
      dE-=site_energy(x_a,y_a,z_a, species_b, CutOff_test);

      dE+=site_energy(x_b,y_b,z_b, species_b, CutOff_test);
      dE-=site_energy(x_b,y_b,z_b, species_a, CutOff_test);
      fprintf(fo,"%d %f \n", CutOff_test, dE);
      dE=0.0;
     }

    fclose(fo);

    // Report on planned move + dE -- for debugging only (makes a ridiculous
    // number of prints...)
    if (DEBUG)
        fprintf(stderr,"MC Move: %d on %d %d %d, %d on %d %d %d --> dE: %f\n",
                lattice[x_a][y_a][z_a], x_a, y_a, z_a,
                lattice[x_b][y_b][z_b], x_b, y_b, z_b,
                dE);

    // OK; now we have dE - proceed with Metropolis Accept/Reject Criteria
    if (dE < 0.0 || exp(-dE * beta) > genrand_real2() )
    {
        //ACCEPT move
        lattice[x_a][y_a][z_a]=species_b; //swap two atoms / species
        lattice[x_b][y_b][z_b]=species_a; 

        ACCEPT++;

        if (EquilibrationChecks) log_dE(dE); 
    }
    else
        //REJECT move; doesn't need to do anything, just update REJECT for
        //stats
        REJECT++;
}
