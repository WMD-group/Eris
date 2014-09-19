/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
void initialise_lattice();

void initialise_lattice()
{
    int x,y,z;
    float angle;

    //Random initial lattice
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                if (genrand_real1()<dipole_fraction) //occupy fraction of sites...
                    random_sphere_point(& lattice[x][y][z]);
                else
                {lattice[x][y][z].x=0.0; lattice[x][y][z].y=0.0; lattice[x][y][z].z=0.0; }
            }
    //Print lattice
    /*
       for (i=0;i<X;i++)
       for (k=0;k<Y;k++)
       printf("\n %f %f %f %f",lattice[i][k].x,lattice[i][k].y,lattice[i][k].z,
       dot(&lattice[i][k],&lattice[i][k]));
       */  
}



