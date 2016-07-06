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
void initialise_lattice_CZTS_supercell();
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



<<<<<<< HEAD
void initialize_lattice_CZTS_supercell()
=======
void initialise_lattice_CZTS_supercell()
>>>>>>> supercell_ordered_initial_lattice
{

// This routine produces a lattice containing a unit cell of CZTS (constructed by inspecting CZTS POSCARs layer by layer)
// Larger systems are then created as supercells of this unit cell using user inputted parameters for X_super, Y_super and Z_super in eris.cfg


// 0 = gap site, 1 = Cu, 2 = Zn, 3 = Sn
// The unit cell is 2x2x4
// Supercells are generated as 2*X_super x 2*Y_super x 4*Z_super

int x,y,z;
    
// Filling corner of the supercell lattice with the unit cell array
<<<<<<< HEAD
  lattice[0][0][0] = 1;
=======
 /* lattice[0][0][0] = 1;
>>>>>>> supercell_ordered_initial_lattice
  lattice[0][0][1] = 0;
  lattice[0][0][2] = 2;
  lattice[0][0][3] = 0;

  lattice[0][1][0] = 0;
  lattice[0][1][1] = 1;
  lattice[0][1][2] = 0;
  lattice[0][1][3] = 3;

  lattice[1][0][0] = 0;
  lattice[1][0][1] = 3;
  lattice[1][0][2] = 0;
  lattice[1][0][3] = 1;

  lattice[1][1][0] = 2;
  lattice[1][1][1] = 0;
  lattice[1][1][2] = 1;
  lattice[1][1][3] = 0;
<<<<<<< HEAD


=======
*/

/*
// Defining unit cell (top left corner of eris terminal display)
// z = 0 layer
lattice[0][0][0] = 2;
lattice[1][0][0] = 0;
lattice[0][1][0] = 0;
lattice[1][1][0] = 1;
// z = 1 layer
lattice[0][0][1] = 0;
lattice[1][0][1] = 1;
lattice[0][1][1] = 3;
lattice[1][1][1] = 0;
// z = 2 layer
lattice[0][0][2] = 1;
lattice[1][0][2] = 0;
lattice[0][1][2] = 0;
lattice[1][1][2] = 2;
// z = 3 layer
lattice[0][0][3] = 0;
lattice[1][0][3] = 3;
lattice[0][1][3] = 1;
lattice[1][1][3] = 0;
*/



// Loops to fill supercell using above unit cell
for (z=0; z<Z_super; z++)
{
    for (y=0; y<Y_super; y++)
    {
        for (x=0; x<X_super; x++)
        {
// z = 0 layer
lattice[0+2*x][0+2*y][3+4*z] = 0;
lattice[1+2*x][0+2*y][3+4*z] = 2;
lattice[0+2*x][1+2*y][3+4*z] = 1;
lattice[1+2*x][1+2*y][3+4*z] = 0;
// z = 1 layer
lattice[0+2*x][0+2*y][2+4*z] = 1;
lattice[1+2*x][0+2*y][2+4*z] = 0;
lattice[0+2*x][1+2*y][2+4*z] = 0;
lattice[1+2*x][1+2*y][2+4*z] = 3;
// z = 2 layer
lattice[0+2*x][0+2*y][1+4*z] = 0;
lattice[1+2*x][0+2*y][1+4*z] = 1;
lattice[0+2*x][1+2*y][1+4*z] = 2;
lattice[1+2*x][1+2*y][1+4*z] = 0;
// z = 3 layer
lattice[0+2*x][0+2*y][0+4*z] = 3;
lattice[1+2*x][0+2*y][0+4*z] = 0;
lattice[0+2*x][1+2*y][0+4*z] = 0;
lattice[1+2*x][1+2*y][0+4*z] = 1;

        }
    }
}




/*
  // Checking initial lattice size
  //int elements = sizeof(lattice)/sizeof(lattice[0][0][0]);
  //int xdim= X_super*2, ydim=Y_super*2, zdim=Z_super*4;
  //fprintf(stderr, "Size of lattice: %d, xdim: %d, ydim: %d, zdim: %d \n", elements, xdim, ydim, zdim);

>>>>>>> supercell_ordered_initial_lattice
  // Creating supercell copies of unit cell in x-direction
  for (x=1;x<X_super;x++)
  {
    lattice[0+2*x][0][0] = lattice[0][0][0];
    lattice[0+2*x][0][1] = lattice[0][0][1];
    lattice[0+2*x][0][2] = lattice[0][0][2];
    lattice[0+2*x][0][3] = lattice[0][0][3];

<<<<<<< HEAD
=======
//debugging
//int xdim1 = 0+2*x; 
//fprintf(stderr, "x supercell1 index: %d \n", xdim1); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[0+2*x][1][0] = lattice[0][1][0];
    lattice[0+2*x][1][1] = lattice[0][1][1];
    lattice[0+2*x][1][2] = lattice[0][1][2];
    lattice[0+2*x][1][3] = lattice[0][1][3];

    lattice[1+2*x][0][0] = lattice[1][0][0];
    lattice[1+2*x][0][1] = lattice[1][0][1];
    lattice[1+2*x][0][2] = lattice[1][0][2];
    lattice[1+2*x][0][3] = lattice[1][0][3];

    lattice[1+2*x][1][0] = lattice[1][1][0];
    lattice[1+2*x][1][1] = lattice[1][1][1];
    lattice[1+2*x][1][2] = lattice[1][1][2];
    lattice[1+2*x][1][3] = lattice[1][1][3];
<<<<<<< HEAD
=======


    //debugging
    //int xdim=1+2*x;
    //fprintf(stderr, "Current xdim: %d, current iteration: %d \n", xdim, x);

>>>>>>> supercell_ordered_initial_lattice
  }

  // Creating supercell copies of unit cell in y-direction
  for (y=1;y<Y_super;y++)
  {
    lattice[0][0+2*y][0] = lattice[0][0][0];
    lattice[0][0+2*y][1] = lattice[0][0][1];
    lattice[0][0+2*y][2] = lattice[0][0][2];
    lattice[0][0+2*y][3] = lattice[0][0][3];

<<<<<<< HEAD
=======
int ydim1 = 0+2*y; 
fprintf(stderr, "y supercell1 index: %d \n", ydim1); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[0][1+2*y][0] = lattice[0][1][0];
    lattice[0][1+2*y][1] = lattice[0][1][1];
    lattice[0][1+2*y][2] = lattice[0][1][2];
    lattice[0][1+2*y][3] = lattice[0][1][3];

<<<<<<< HEAD
=======

int ydim1_again = 0+2*y; 
fprintf(stderr, "y supercell1 index again: %d \n", ydim1_again); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[1][0+2*y][0] = lattice[1][0][0];
    lattice[1][0+2*y][1] = lattice[1][0][1];
    lattice[1][0+2*y][2] = lattice[1][0][2];
    lattice[1][0+2*y][3] = lattice[1][0][3];

<<<<<<< HEAD
=======
int ydim2 = 1+2*y; 
fprintf(stderr, "y supercell2 index: %d \n", ydim2); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[1][1+2*y][0] = lattice[1][1][0];
    lattice[1][1+2*y][1] = lattice[1][1][1];
    lattice[1][1+2*y][2] = lattice[1][1][2];
    lattice[1][1+2*y][3] = lattice[1][1][3];
<<<<<<< HEAD
  }
=======


int ydim2_again = 1+2*y; 
fprintf(stderr, "y supercell2 index again: %d \n", ydim2_again); 

    //debugging
    //int ydim=1+2*y;
    //fprintf(stderr, "Current ydim: %d, current iteration: %d \n", ydim, y);
    
  }

>>>>>>> supercell_ordered_initial_lattice
  // Creating supercell copies of unit cell in z-direction
  for (z=1;z<Z_super;z++)
  {
    lattice[0][0][0+4*z] = lattice[0][0][0];
    lattice[0][0][1+4*z] = lattice[0][0][1];
    lattice[0][0][2+4*z] = lattice[0][0][2];
    lattice[0][0][3+4*z] = lattice[0][0][3];

<<<<<<< HEAD
=======
int zdim1 = 3+4*z; 
fprintf(stderr, "z supercell1 index: %d \n", zdim1); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[0][1][0+4*z] = lattice[0][1][0];
    lattice[0][1][1+4*z] = lattice[0][1][1];
    lattice[0][1][2+4*z] = lattice[0][1][2];
    lattice[0][1][3+4*z] = lattice[0][1][3];

<<<<<<< HEAD
=======

int zdim1_again = 3+4*z; 
fprintf(stderr, "z supercell1 index again: %d \n", zdim1_again); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[1][0][0+4*z] = lattice[1][0][0];
    lattice[1][0][1+4*z] = lattice[1][0][1];
    lattice[1][0][2+4*z] = lattice[1][0][2];
    lattice[1][0][3+4*z] = lattice[1][0][3];

<<<<<<< HEAD
=======

int zdim2 = 3+4*z; 
fprintf(stderr, "z supercell2 index: %d \n", zdim2); 

>>>>>>> supercell_ordered_initial_lattice
    lattice[1][1][0+4*z] = lattice[1][1][0];
    lattice[1][1][1+4*z] = lattice[1][1][1];
    lattice[1][1][2+4*z] = lattice[1][1][2];
    lattice[1][1][3+4*z] = lattice[1][1][3];
<<<<<<< HEAD
  }
=======


int zdim2_again = 3+4*z; 
fprintf(stderr, "z supercell2 index again: %d \n", zdim2_again); 


// Printing full array for debugging
for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
                   fprintf(stderr,"%d ", lattice[i][j][0]);
        }
        fprintf(stderr,"\n");
}

    //debugging
    //int zdim=1+4*z;
    //fprintf(stderr, "Current zdim: %d, current iteration: %d \n", zdim, z);
  
  }
*/


>>>>>>> supercell_ordered_initial_lattice

}


void initialise_lattice_CZTS_randomized()
{
    initialise_lattice_CZTS_supercell(); // start with correct stoichometry

    // Randomizing CZTS lattice by swapping species
    int shuffles=10; // Defining number of swap attempts based on tota no. of sites
    int i;
    int x_swap, y_swap, z_swap;

    int x,y,z;
    int species;

    for (i=0;i<shuffles;i++)
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                if (lattice[x][y][z]==0) continue;
                species=lattice[x][y][z];

                // Choosing two random sites to swap
                x_swap = rand_int(X);
                y_swap = rand_int(Y);
                z_swap = rand_int(Z);

                // Checking that neither species is a gap (do not want to randomize gap sites)
                if (lattice[x_swap][y_swap][z_swap]==0) continue;

                if (freezeSn) // don't swap Tin if requested
                    if (species==Sn || lattice[x_swap][y_swap][z_swap]==Sn)
                        continue;
      
                // Swapping species on site1 with species on site2 and vice versa (using temporary intermediate variable)
                lattice[x][y][z]=lattice[x_swap][y_swap][z_swap];
                lattice[x_swap][y_swap][z_swap]=species;
            }

}
