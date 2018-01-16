/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
static double potential_at_site(int x, int y, int z);
static double potential_at_site_r_test(int x, int y, int z, int r_cutoff);
static void lattice_potential_log(FILE *log);
void lattice_potential_XY(char * filename);
void lattice_potential_XYZ(char * filename);
void equil_lattice_potential(char * filename);
void lattice_potential_r_test(char * filename);
static double lattice_energy_log(FILE *log);
double landau_order();

void outputpotential_png(char * filename);
void radial_distribution_function(char * filename, int speciesA, int speciesB);
void radial_distribution_function_allsites();
void radial_distribution_function_allsites_initial(); 

void outputlattice_xyz(char * filename);
void outputlattice_xyz(char * filename);
void outputlattice_xyz_overprint(char * filename);
void outputlattice_pymol_cgo(char * filename);
void outputlattice_dumb_terminal();

void outputlattice_stoichometry();

void generate_gulp_input(int temp, char * filename);  
static void log_dE(float dE);



static double potential_at_site(int x, int y, int z) 
{
    int dx,dy,dz=0;
    double pot=0.0;
    float d;
    struct dipole r;

    for (dx=-POTENTIAL_CUTOFF;dx<POTENTIAL_CUTOFF;dx++)
        for (dy=-POTENTIAL_CUTOFF;dy<POTENTIAL_CUTOFF;dy++)
            for (dz=-POTENTIAL_CUTOFF;dz<POTENTIAL_CUTOFF;dz++)
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                r.x=(float)(dx); r.y=(float)(dy); r.z=(float)(dz);

                d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

                if (d>(float)POTENTIAL_CUTOFF) continue; // Cutoff in d

                // pot(r) = 1/4PiEpsilon * p.r / r^3
               
               //for CZTS
                int species;
                species=lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];
                double q;
                q=FormalCharge[species];

                pot+= q/(double)d;
            }
    return(pot);
}


static double potential_at_site_r_test(int x, int y, int z, int r_cutoff) 
{
    int dx,dy,dz=0;
    double pot=0.0;
    float d;
    struct dipole r;

    for (dx=-r_cutoff;dx<r_cutoff;dx++)
        for (dy=-r_cutoff;dy<r_cutoff;dy++)
            for (dz=-r_cutoff;dz<r_cutoff;dz++)
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                r.x=(float)(dx); r.y=(float)(dy); r.z=(float)(dz);

                d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

                if (d>(float)r_cutoff) continue; // Cutoff in d

                // pot(r) = 1/4PiEpsilon * p.r / r^3
               
               //for CZTS
                int species;
                species=lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z];
                double q;
               // Use effetive charges (defined in eris.cfg) to account for S anion in between each cation
                q=EffectiveCharge[species];

                pot+= q/(double)d;
            }
    return(pot);
}


//Calculates dipole potential along trace of lattice
static void lattice_potential_log(FILE *log)
{
    int x,y,z;
    double pot;

    y=Y/2; //trace across centre of material. I know, I know, PBCs.
    z=0;
    for (x=0;x<X;x++)
    {
        pot=0.0;
        for (y=0;y<Y;y++)
            pot+=potential_at_site(x,y,z);
        fprintf(log,"%d %f %f\n",x,pot/(double)Y,potential_at_site(x,Y/2,z));
    }

}

//Calculates dipole potential across XY lattice
void lattice_potential_XY(char * filename)
{
    int x,y;
    double pot;
    FILE *fo;
    fo=fopen(filename,"w");

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            fprintf(fo,"%d %d %f\n",x,y,potential_at_site(x,y,0));
}

//Calculates dipole potential across XYZ volume
void lattice_potential_XYZ(char * filename)
{
    int x,y,z;
    int atoms;
    double pot,mean,variance;
    FILE *fo;
    fo=fopen(filename,"a");

    FILE *fvariance;
    fvariance=fopen("variance.dat","a"); // append to variance file

    mean=0.0; atoms=0;
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                // log potential at all sites (Cu,Zn,Sn)
                pot=potential_at_site(x,y,z);
                fprintf(fo,"%d %d %d %d %f\n",lattice[x][y][z],x,y,z,pot);
                
                if (lattice[x][y][z]==3) // only count tin towards mean / variance
                {
                    mean+=pot;
                    atoms++;
                }
            }
    mean/=atoms;

    variance=0.0; atoms=0;
     for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                if (lattice[x][y][z]==3)
                {
                    pot=potential_at_site(x,y,z);
                    variance+=(pot-mean)*(pot-mean);
                    atoms++;
                }
            }
    variance/=atoms;

    fprintf(stderr,"T: %04d Mean: %f Variance(rigorous): %f TinAtoms: %d Total(X*Y*Z):%d\n",
            T,mean,variance,atoms,X*Y*Z);
    fprintf(fo,"# T: %04d Mean: %f Variance(rigorous): %f TinAtoms: %d Total(X*Y*Z):%d\n",
            T,mean,variance,atoms,X*Y*Z);
    fclose(fo);
    
    fprintf(fvariance,"T: %04d Mean: %f Variance(rigorous): %f TinAtoms: %d Total(X*Y*Z):%d\n",
            T,mean,variance,atoms,X*Y*Z);
    fclose(fvariance);
}

// Outputs PNG (bitmap picture) for the potential across the X-Y plane at Z=0
void outputpotential_png(char * filename)
{
    int i,k,pixel;
    FILE *fo;
    fo=fopen(filename,"w");

    fprintf (fo,"P2\n%d %d\n%d\n", X, Y, SHRT_MAX);

    for (i=0;i<X;i++)
    {
        for (k=0;k<Y;k++)
        {
            pixel=SHRT_MAX/2+(int)(SHRT_MAX*0.1*potential_at_site(i,k,0));

            // Bounds checking :^)
            if (pixel<0) pixel=0;
            if (pixel>SHRT_MAX) pixel=SHRT_MAX; // PNG P2 bitmap has binary values from 0..SHRT_MAX

            fprintf(fo,"%d ",pixel);
        }
        fprintf(fo,"\n");
    }

}

// Generates RDFs for {All Species} <-> {All Species}, constructing the
// filename then calling the radial_distribution_function which does the actual
// work
void radial_distribution_function_allsites()
{
    int speciesA,speciesB;
// Variable to contain filename for RDF output
    char RDF_filename[100]; 
    
    for (speciesA=1;speciesA<=3;speciesA++)
        for (speciesB=speciesA;speciesB<=3;speciesB++)
        {
            sprintf(RDF_filename,"RDF_%c_%c_Temp_%04d.dat",specieslookup[speciesA],specieslookup[speciesB],T); // automatically construct filename
            radial_distribution_function(RDF_filename, speciesA, speciesB); // run the RDF for this pair

            fprintf(stderr,"RDF: speciesA: %d speciesB: %d filename: %s\n",speciesA,speciesB,RDF_filename); // for debugging
        }
}           

// Extra wrapper function added for outputting initial RDF to a separate file before performing any MC steps
void radial_distribution_function_allsites_initial() 
{
    int speciesA,speciesB;
// Variable to contain filename for RDF output
    char RDF_filename[100]; 
    
    for (speciesA=1;speciesA<=3;speciesA++)
        for (speciesB=speciesA;speciesB<=3;speciesB++)
        {
            sprintf(RDF_filename,"RDF_%c_%c_initial.dat",specieslookup[speciesA],specieslookup[speciesB]); // automatically construct filename 
            radial_distribution_function(RDF_filename, speciesA, speciesB); // run the RDF for this pair

            fprintf(stderr,"RDF: speciesA: %d speciesB: %d filename: %s\n",speciesA,speciesB,RDF_filename); // for debugging
        }
}

void radial_distribution_function(char * filename, int speciesA, int speciesB )
// Calculates RDF for on-lattice material
//    arguments are filename to write to, then atomic-pair to calculate RDF for
{
    int x,y,z;
    int dx,dy,dz;
    int i;
    FILE *fo;
    fo=fopen(filename, "w"); // w setting used so RDF file is overwritten, not appended

    int distance_squared;

    const int CUTOFF=9;

    struct dipole n;
    float d;
    float pair_correlation;

    // define data structures to keep histogram counts in
    float RDF[(CUTOFF*CUTOFF)+1];
    int site_count[(CUTOFF*CUTOFF)+1];
    for (i=0;i<=CUTOFF*CUTOFF;i++) // Zero histogram arrays
        { RDF[i]=0.0; site_count[i]=0; }

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                for (dx=-CUTOFF;dx<=CUTOFF;dx++)
                    for (dy=-CUTOFF;dy<=CUTOFF;dy++)
                       for (dz=-CUTOFF;dz<=CUTOFF;dz++)
                        {
                            distance_squared=dx*dx + dy*dy + dz*dz;
                            if (distance_squared>CUTOFF*CUTOFF) continue; // skip ones that exceed spherical limit of CUTOFF

                            if (lattice[x][y][z]!=speciesA) continue; // if not speciesA, discard this data

                            pair_correlation =  lattice[(x+dx+X)%X][(y+dy+Y)%Y][(z+dz+Z)%Z]==speciesB ? 1.0 : 0.0; 
                            //complicated modulus arithmatic deals with PBCs

                            // Collect sum into RDF histogram 
                            RDF[distance_squared]+=pair_correlation;
                            site_count[distance_squared]++; // count of number of lattice sites at this separation, for normalisation purposes
                        }

    // Weight counts into a RDF & output to file called RDF_T_temp.dat as specified in eris-main.c
    fprintf(fo,"# r^2 r RDF(lattice) RDF(r^2) \t RDF[notnormalised] SitesInspected[r^2] T\n");
    for (i=0;i<CUTOFF*CUTOFF;i++)
    {
        if (site_count[i]>0) //i.e. any data here...
        {
            fprintf(fo, "%d %f %f %f \t %d \t %d \t %d\n",i,sqrt(i),
                    RDF[i]/site_count[i],
                    RDF[i]/(float)i,
                    (int)RDF[i],site_count[i],T);
        }
    }
    fprintf(fo,"\n"); //starts as new dataset in GNUPLOT --> discontinuous lines

    fclose(fo);
}

void outputlattice_xyz(char * filename)
{
    int i,j,k;
    FILE *fo;
    const char * atom[] = {
            "Nu",
            "Cu",
            "Zn",
            "Sn",
            "Nu"
    };
    const float d=2.72; // Angstrom spacing of lattice to map to real space coords

    fo=fopen(filename,"w");
   
    fprintf(fo,"%d\n\n",X*Y*Z);

    for (i=0;i<X;i++)
        for (j=0;j<Y;j++)
            for (k=0;k<Z;k++)
                fprintf(fo,"%s %f %f %f\n",atom[lattice[i][j][k]],d*(float)i,d*(float)j,d*(float)k);
    fclose(fo);
}

#define ZSCALE 5.0 // Scales Z-axis in Pymol xyz / CGO outputs

float DMAX=55.0; //sensible starting value...
float DMEAN=0.0;

// Prints 2D reproresentation of the lar_cutof using ANSI
// colours to make it look pretty; and presenting it side-by-side with the
// potential map for the simulated volume
void outputlattice_dumb_terminal()
{
    int x,y;
    float a;
    int z=0;
    float new_DMAX=0.0; //used to calibrate next colour scale, based on present maxima of data
    float variance=0.0; // sum of potential^2
    float mean=0.0;
    float potential;

    fprintf(stderr,"%*s%*s\n",X+3, "SPECIES", (2*X)+4,"POTENTIAL"); //padded labels


    for (z=0;z<DumbTerminalLayers;z++) // number of layers to display on Z
    {
        fprintf(stderr,"Z=%d\n",z);

        for (y=0;y<Y;y++)
            for (x=0;x<X;x++)
            {
                potential=potential_at_site(x,y,z);
                if (fabs(potential-DMEAN)>new_DMAX)
                    new_DMAX=fabs(potential-DMEAN); // used to calibrate scale - technically this changes
            }
        DMAX=new_DMAX;

        for (y=0;y<Y;y++)
        {
            for (x=0;x<X;x++)
            {
                a=lattice[x][y][z];

                fprintf (stderr,"%c[%d",27,31+((int)a)%8 ); // Sets colour of output routine
                fprintf(stderr,";7"); // inverted colours
                char species=specieslookup[(int)a];

                fprintf(stderr,"m%c %c[0m",species,27);  // prints 1-character reference to species
                fprintf(stderr,"%c[37m%c[0m",27,27); //RESET
            }

            // OK - now potential plot :^)
            //        const char * density=".,:;o*O#"; //increasing potential density
            const char * density="012345689";
            fprintf(stderr,"    ");
            for (x=0;x<X;x++)
            {
                potential=potential_at_site(x,y,z);

                variance+=potential*potential;
                mean+=potential;

                if (fabs(potential-DMEAN)>new_DMAX)
                    new_DMAX=fabs(potential-DMEAN); // used to calibrate scale - technically this changes
                //printf("%f\t",potential); //debug routine to get scale

                //fprintf(stderr,"%c[%d",27,31+((int)(8.0*fabs(potential)/DMAX))%8); //8 colours
                //fprintf(stderr,"%c[48;5;%d",27,17+(int)(214.0*fabs(potential)/DMAX)); // Xterm 256 color map - (16..231)
                fprintf(stderr,"%c[48;5;%d",27,232+12+(int)(12.0*(potential-DMEAN)/DMAX)); // Xterm 256 color map - shades of grey (232..255)
                // https://code.google.com/p/conemu-maximus5/wiki/AnsiEscapeCodes#xterm_256_color_processing_requirements

                //if (potential<0.0) // if negative
                //    fprintf(stderr,";7"); // bold

                a=lattice[x][y][z];

                char species=specieslookup[(int)a];  // prints 1-character reference to species

                fprintf(stderr,"m%c%c%c[0m",density[(int)(8.0*fabs(potential-DMEAN)/DMAX)],species,27);
            }

            fprintf(stderr,"\n");
        }
        mean=mean/(X*Y);
        DMEAN=mean; // for calibration of scale

        variance=variance/(X*Y); 
        fprintf(stdout,"T: %d DMAX: %f new_DMAX: %f (not quite) variance: %f mean: %f\n",T,DMAX,new_DMAX,variance,mean);
    }
    DMAX=new_DMAX; // infinite fast following - but leads to fluctuations at steady state
    if (DMAX==0.0) DMAX=1.0; //avoid divide by zero for all-zero pot
}

void lattice_energy ()
{
    int dx,dy,dz;

    int x=5,y=5,z=4; //which site to move around

    int DELTA=3;;
    for (dx=0;dx<DELTA;dx++)
        for (dy=0;dy<DELTA;dy++)
            for (dz=0;dz<DELTA;dz++)
            {
                if (lattice[(x+dx+X)%X][(y+dy+X)%X][(z+dz+Z)%Z]==0)
                    continue;

                printf("X:%d Y:%d Z:%d dx: %d dy: %d dz: %d",x,y,z,dx,dy,dz);
                //int cutoff=3;
                for (int cutoff=1;cutoff<=10;cutoff++)
                {
                    double site_E=site_energy(x,y,z,lattice[x][y][z],cutoff);
                    printf("\n\t r: %d E:%f",cutoff,site_E); 
                    
                    double dE=0.0; 
                    int species_a=lattice[x][y][z]; 
                    int species_b=lattice[(x+dx+X)%X][(y+dy+X)%X][(z+dz+Z)%Z];
                    
//                    dE+=site_energy(x,y,z, species_a, cutoff);
//                    dE-=site_energy(x,y,z, species_b, cutoff);
//                    dE+=site_energy((x+dx+X)%X,(y+dy+Y)%Y,(z+dz+Z)%Z, species_b, cutoff);
//                    dE-=site_energy((x+dx+X)%X,(y+dy+Y)%Y,(z+dz+Z)%Z, species_a, cutoff);
 
                    dE+=site_energy_stencil(x,y,z, species_a, cutoff, x,y,z);
                    dE-=site_energy_stencil(x,y,z, species_b, cutoff, x,y,z);
                    dE+=site_energy_stencil((x+dx+X)%X,(y+dy+Y)%Y,(z+dz+Z)%Z, species_b, cutoff, x, y, z);
                    dE-=site_energy_stencil((x+dx+X)%X,(y+dy+Y)%Y,(z+dz+Z)%Z, species_a, cutoff, x, y, z);
                    
                    printf("\tdE:%f",dE); 
                }
                printf("\n");
            }
}

void outputlattice_stoichometry()
{
    int x,y,z,i,a;
   
    // data structure for histogram of species count
    int histogram[5];
    for (i=0;i<5;i++)
        histogram[i]=0;

    // loop over lattice and sum histogram
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            histogram[lattice[x][y][z]]++;

    // output histogram sums for user to read
    printf("Lattice stoichometry:\n");
    for (i=0;i<5;i++)
        printf("    Species: %c Counts: %d\n",specieslookup[i],histogram[i]);

    if (histogram[1]!=2*histogram[2] || histogram[1]!=2*histogram[3])
        fprintf(stderr,"Warning! Lattice not stoichometric!\n");

/* Added by Suzy to immediately exit if lattice not stoichometric   
    // conditional statements to check for stoichiometry and exit with error message if off-stoichiometry
    if (histogram[0] != 2*histogram[1])
    {
        printf("Lattice not stoichiometric (number of gaps compared to Cu and maybe more problems!) \n");
        exit(EXIT_FAILURE);
    }
 
    if (histogram[1] != 2*histogram[2])
    {
        printf("Lattice not stoichiometric (number of Cu compared to Zn and maybe compared to Sn!) \n");
        exit(EXIT_FAILURE);   
    }
    if (histogram[1] != 2*histogram[3])
    {
        printf("Lattice not stoichiometric (number of Cu compared to Sn) \n");
        exit(EXIT_FAILURE);    
    }
    //printf("%d \n",histogram[1]);
*/
}

void equil_lattice_potential(char * filename)
{
    int x,y,z;
    double pot;

    FILE *fo;
    fo=fopen(filename,"w"); 

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                // log potential at only Sn sites
                if (lattice[x][y][z]==3)
                {
                  pot=potential_at_site(x,y,z);
                  fprintf(fo, "%f \n",pot);
                }                

            }
fclose(fo);
}

void lattice_potential_r_test(char * filename)
{
    int x,y,z,r_cutoff;
    double pot;

    // Defining max r_cutoff for potential calculation based on lattice dimensions
    int CutOffMax = floor(min(X,Y,Z)/2);

    FILE *fo;
    fo=fopen(filename,"w"); 

    for (r_cutoff=1;r_cutoff<CutOffMax;r_cutoff++)
    {
      for (x=0;x<X;x++)
          for (y=0;y<Y;y++)
              for (z=0;z<Z;z++)
              {
                  // log potential at only Sn sites
                  if (lattice[x][y][z]==3)
                  {
                    pot=potential_at_site_r_test(x,y,z,r_cutoff);
                    fprintf(fo, "%f \n",pot);
                  }                

              }
    }  
fclose(fo);
}

void generate_gulp_input(int temp, char * filename)
{
    int i,j,k;
    char selected_site[100];
    FILE *fo;
    const char * atom[] = {
            "Nu",
            "Cu",
            "Zn",
            "Sn",
            "Nu"
    };
    
    const char * formal_charge[] = {
            "empty site",
            "1.0",
            "2.0",
            "4.0",
            "empty site again"
    };
    const float d=2.72; // Angstrom spacing of lattice to map to real space coords

    fo=fopen(filename,"w");
   
//    fprintf(fo,"%d\n\n",X*Y*Z);

    // Writing top lines of gulp input file
    const float X_dim=d*X, Y_dim=d*Y, Z_dim=d*Z;

    fprintf(fo, "# Keywords: \n");
    fprintf(fo, "# \n");
    fprintf(fo, "pot \n");
    fprintf(fo, "# \n");
    fprintf(fo, "# Options: \n");
    fprintf(fo, "# \n");
    fprintf(fo, "cell \n");
    fprintf(fo, "%f %f %f 90.000000 90.000000 90.000000 \n", X_dim, Y_dim, Z_dim);
    fprintf(fo, "cartesian \n");

        // Adding S anions to top of the coordinates list based on the fixed S positions in a unit cell, expanded using supercell parameters
        for (i=0; i<X/2; i++)
        {
          for (j=0; j<Y/2; j++)
          {
            for (k=0; k<Z/4; k++)
            {
              // x- and y- coordinates are swapped to make visual of xyz file consistent with POSCAR in VESTA
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.758700013*2.0*d)+(2.0*d*i), (0.254750013*2.0*d)+(2.0*d*j), (0.877870023*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.241300002*2.0*d)+(2.0*d*i), (0.745249987*2.0*d)+(2.0*d*j), (0.877870023*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.254750013*2.0*d)+(2.0*d*i), (0.241300002*2.0*d)+(2.0*d*j), (0.122129999*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.745249987*2.0*d)+(2.0*d*i), (0.758700013*2.0*d)+(2.0*d*j), (0.122129999*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.258700013*2.0*d)+(2.0*d*i), (0.754750013*2.0*d)+(2.0*d*j), (0.377869993*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.741299987*2.0*d)+(2.0*d*i), (0.245249987*2.0*d)+(2.0*d*j), (0.377869993*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.754750013*2.0*d)+(2.0*d*i), (0.741299987*2.0*d)+(2.0*d*j), (0.622129977*4.0*d)+(4.0*d*k));
              fprintf(fo,"S core %f %f %f -2.0 \n", (0.245249987*2.0*d)+(2.0*d*i), (0.258700013*2.0*d)+(2.0*d*j), (0.622129977*4.0*d)+(4.0*d*k));
              
            }
          }
        }


        // Looping over lattice to write coordinates of cations (C, Z, T) to gulp input file
        for (i=0;i<X;i++)
          for (j=0;j<Y;j++)
            for (k=0;k<Z;k++)

                if (lattice[i][j][k]==0) continue; //avoid writing gap sites to gulp input file
                else fprintf(fo,"%s core %f %f %f %s \n",atom[lattice[i][j][k]],d*(float)i,d*(float)j,d*(float)k,formal_charge[lattice[i][j][k]]);
   
    fprintf(fo, "output xyz final_lattice_T_%04d.xyz",T);
   // fprintf(fo, "space \n");
    // fprintf(fo, "82 \n");
     fclose(fo);

}

// dE log to compare to Boltzmann stats below here 
#define dE_BINS 100
int endo_bins[dE_BINS], exo_bins[dE_BINS];

double sum_dE=0.0;

static void reset_dE()
{
    sum_dE=0.0;
    for (int i=0;i<dE_BINS;i++)
        endo_bins[i]=exo_bins[i]=0;
}

static void report_dE()
{
         printf("Change in energy: %f\n",sum_dE);

         printf("Histogram of exo/endo: \n");
         for (int i=0; i<dE_BINS; i++)
         {
             printf("%d",endo_bins[i]);
             printf("/%d ",exo_bins[i]);
             endo_bins[i]=0;
             exo_bins[i]=0;
         }
         printf("\n");
}

static void log_dE(float dE)
{
     sum_dE+=(double)dE;
     int my_bin=(int) (pow(fabsf(dE),0.5)/0.05); // sqrt(data) for storage

     if (my_bin>=dE_BINS) my_bin=dE_BINS-1; // range checking logic to avoid segfaults

//     printf("beta: %f dE: %f logf(fabsf(dE)): %f logf(fabsf(dE*beta)): %f\n",beta,dE,logf(fabsf(dE)),pow(fabsf(dE),0.5));
     if (dE>0.0)
        endo_bins[ my_bin ] ++;
     else
         exo_bins[ my_bin ] ++;

}

