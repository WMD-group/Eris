/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
void initialise_lattice_random();
void initialise_lattice_stripe();
void initialise_lattice_CZTS();
void initialise_lattice_CZTS_randomized();

void initialise_lattice_random()
{
    int x,y,z;
    int species;

    //Random initial lattice
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                species=rand_int(4)+1;
                // 0: Nothing
                // 1: Copper (I)
                // 2: Zinc  (II)
                // 3: Tin   (IIII)
                // 4: Copper (Dummy --> change to 4)
                if (species==4) species=1; //Twice as much copper
            
                if ((x+y+z)%2==0)
                    lattice[x][y][z]=species; 
                else
                    lattice[x][y][z]=0; // gaps for the FCC sublattice
            }
    //Print lattice
    /*
       for (i=0;i<X;i++)
       for (k=0;k<Y;k++)
       printf("\n %f %f %f %f",lattice[i][k].x,lattice[i][k].y,lattice[i][k].z,
       dot(&lattice[i][k],&lattice[i][k]));
       */  

}

void initialise_lattice_stripe()
{
    int x,y,z;
    int species;

    //Striped lattice...
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                species=1+ ((y*4)/Y);
                // 0: Nothing
                // 1: Copper (I)
                // 2: Zinc  (II)
                // 3: Tin   (IIII)
                // 4: Copper (Dummy --> change to 4)
               if (species==3) species=1; //Twice as much copper

                lattice[x][y][z]=species; 
            }
}

void initialise_lattice_CZTS()
{
    int x,y,z;
    int species;

    //Striped lattice...
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                species=1+ (x+2*y+2*z)%4;

                // 0: Nothing
                // 1: Copper (I)
                // 2: Zinc  (II)
                // 3: Tin   (IIII)
                // 4: Copper (Dummy --> change to 4)
                if (species==4) species=1; //Twice as much copper

                if ((x+y+z)%2==0)
                    lattice[x][y][z]=species; 
                else
                    lattice[x][y][z]=0; // gaps for the FCC sublattice
            }
 
}




void initialise_lattice_CZTS_randomized()
{
    int x,y,z;
    int species;

    initialise_lattice_CZTS(); // start with correct stoichometry

    // Randomizing CZTS lattice by swapping species
    int swaps=X*Y*Z*1000; // Defining number of swap attempts based on tota no. of sites
    int i;
    int x_swap1, y_swap1, z_swap1, x_swap2, y_swap2, z_swap2, tmp_species;
    for (i=0;i<swaps;i++)
    {
      // Choosing two random sites to swap
      x_swap1 = rand_int(X);
      y_swap1 = rand_int(Y);
      z_swap1 = rand_int(Z);

      x_swap2 = rand_int(X);
      y_swap2 = rand_int(Y);
      z_swap2 = rand_int(Z);
      // Checking that neither species is a gap (do not want to randomize gap sites)
      if (lattice[x_swap1][y_swap1][z_swap1]==0){
          continue;
      }
      if (lattice[x_swap2][y_swap2][z_swap2]==0){
          continue;
      }
      // Swapping species on site1 with species on site2 and vice versa (using temporary intermediate variable)
      tmp_species = lattice[x_swap1][y_swap1][z_swap1];
      lattice[x_swap1][y_swap1][z_swap1] = lattice[x_swap2][y_swap2][z_swap2];
     lattice[x_swap2][y_swap2][z_swap2] = tmp_species;

    }

}
