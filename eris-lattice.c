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
		// 3: Also Copper... (II)
                // 4: Tin   (IIII)
                if (species==3) species=1; //Twice as much copper

                lattice[x][y][z]=species; 
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
		// 3: Copper (Equiv to 1)
                // 3: Tin   (IIII)
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
		// 3: Copper (equiv to 1)
                // 4: Tin   (IIII)
                if (species==3) species=1; //Twice as much copper

                lattice[x][y][z]=species; 
            }
 
}
