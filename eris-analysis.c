/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

// Prototypes...
static void lattice_angle_log(FILE *log);
static double polarisation();
static double dipole_potential(int x, int y, int z);
static void lattice_potential_log(FILE *log);
void lattice_potential_XY(char * filename);
void lattice_potential_XYZ(char * filename);
static double lattice_energy_log(FILE *log);
double landau_order();

void outputpotential_png(char * filename);
void radial_distribution_function();

void outputlattice_xyz(char * filename);
void outputlattice_pnm(char * filename);
void outputlattice_ppm_hsv(char * filename);
void outputlattice_svg(char * filename);
void outputlattice_xyz(char * filename);
void outputlattice_xyz_overprint(char * filename);
void outputlattice_pymol_cgo(char * filename);
void outputlattice_dumb_terminal();


//Calculate dipole potential at specific location
static double dipole_potential(int x, int y, int z) 
{
    int dx,dy,dz=0;
    double pot=0.0;
    float d;
    struct dipole r;

    for (dx=-POTENTIAL_CUTOFF;dx<POTENTIAL_CUTOFF;dx++)
        for (dy=-POTENTIAL_CUTOFF;dy<POTENTIAL_CUTOFF;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-POTENTIAL_CUTOFF;dz<POTENTIAL_CUTOFF;dz++)
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                r.x=(float)(dx); r.y=(float)(dy); r.z=(float)(dz);

                d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

                if (d>(float)POTENTIAL_CUTOFF) continue; // Cutoff in d

                 // pot(r) = 1/4PiEpsilon * p.r / r^3
               // pot(r) = 1.0d0/4PiEpsilon * p.r / r^3
                // Electric dipole potential
// FIXME
                 //dipole                pot+=dot(& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z] ,& r)/(d*d*d);
               //for CZTS
                pot+= (float)(lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z]-2)/d;
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
            pot+=dipole_potential(x,y,z);
        fprintf(log,"%d %f %f\n",x,pot/(double)Y,dipole_potential(x,Y/2,z));
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
            fprintf(fo,"%d %d %f\n",x,y,dipole_potential(x,y,0));
}

//Calculates dipole potential across XYZ volume
void lattice_potential_XYZ(char * filename)
{
    int x,y,z;
    int atoms;
    double pot,mean,variance;
    FILE *fo;
    fo=fopen(filename,"w");

    FILE *fvariance;
    fvariance=fopen("variance.dat","a"); // append to variance file

    mean=0.0; atoms=0;
    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                // log potential at all sites (Cu,Zn,Sn)
                pot=dipole_potential(x,y,z);
                fprintf(fo,"%d %d %d %f\n",x,y,z,pot);
                
                if (lattice[x][y][z]==4) // only count tin towards mean / variance
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
                if (lattice[x][y][z]==4)
                {
                    pot=dipole_potential(x,y,z);
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
            pixel=SHRT_MAX/2+(int)(SHRT_MAX*0.1*dipole_potential(i,k,0));

            // Bounds checking :^)
            if (pixel<0) pixel=0;
            if (pixel>SHRT_MAX) pixel=SHRT_MAX;

            fprintf(fo,"%d ",pixel);
        }
        fprintf(fo,"\n");
    }

}

void radial_distribution_function()
// Calculates RDF for on-lattice material
// Currently prints to stdout
{
    int x,y,z;
    int dx,dy,dz;
    int i;

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

                            if (lattice[x][y][z]!=2) continue; // if not Copper

                            //complicated modulus arithmatic deals with PBCs
                            pair_correlation = lattice[x][y][z]==lattice[(x+dx+X)%X][(y+dy+Y)%Y][(z+dz+Z)%Z] ? 1.0 : 0.0;
                            // if species at site A = species at site B; add one to the pair correlation
                            // correlation...

                            // Collect sum into RDF histogram 
                            RDF[distance_squared]+=pair_correlation;
                            site_count[distance_squared]++; // count of number of lattice sites at this separation, for normalisation purposes
                        }

    // Weight counts into a RDF & output to STDOUT
    printf("# r^2 r RDF(lattice) RDF(r^2) \t RDF[notnormalised] SitesInspected[r^2] T\n");
    for (i=0;i<CUTOFF*CUTOFF;i++)
    {
        if (site_count[i]>0) //i.e. any data here...
        {
            printf("%d %f %f %f \t %d \t %d \t %d\n",i,sqrt(i),
                    RDF[i]/site_count[i],
                    RDF[i]/(float)i,
                    (int)RDF[i],site_count[i],T);
        }
    }
    printf("\n"); //starts as new dataset in GNUPLOT --> discontinuous lines

    return; //Ummm
}

void outputlattice_xyz(char * filename)
{
    int i,j,k;
    FILE *fo;
    const char * atom[] = {
            "Nu",
            "Cu",
            "Zn",
            "Nu",
            "Sn"
    };
    const float d=3.8; // Angstrom spacing of lattice to map to real space coords

    fo=fopen(filename,"w");
    
    for (i=0;i<X;i++)
        for (j=0;j<Y;j++)
            for (k=0;k<Z;k++)
                fprintf(fo,"%s %f %f %f\n",atom[lattice[i][j][k]],d*(float)i,d*(float)j,d*(float)k);
    fclose(fo);
}

void outputlattice_png(char * filename) // TODO: Fix me for new code
{
    int i,k;
    FILE *fo;
    fo=fopen(filename,"w");

    fprintf (fo,"P2\n%d %d\n%d\n", X, Y, SHRT_MAX);

    for (i=0;i<X;i++)
    {
        for (k=0;k<Y;k++)
//            fprintf(fo,"%d ",(int)(SHRT_MAX*atan2(lattice[i][k][0].y,lattice[i][k][0].x)/(2*M_PI)));
        fprintf(fo,"\n");
    }

}

// Outputs a PPM bitmap of lattice dipole orientation on a HSV colourwheel
void outputlattice_ppm_hsv(char * filename)
{
    int i,k;
    float angle;

    float r,g,b; // RGB
    float h,s,v; // HSV
    float p,t,q,f; // intemediates for HSV->RGB conversion
    int hp;

    FILE *fo;
    fo=fopen(filename,"w");

    //Set Saturation + Value, vary hue
    s=0.6; v=0.8;

    fprintf (fo,"P6\n%d %d\n255\n", X, Y);

    for (i=0;i<X;i++) //force same ordering as SVG...
        for (k=0;k<Y;k++)
        {
/* TODO: Fix for new lattice types
            h=M_PI+atan2(lattice[i][k][0].y,lattice[i][k][0].x); //Nb: assumes 0->2PI interval!
            v=0.5+0.4*lattice[i][k][0].z; //darken towards the south (-z) pole
            s=0.6-0.6*fabs(lattice[i][k][0].z); //desaturate towards the poles
*/
            // http://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
            hp=(int)floor(h/(M_PI/3.0)); //radians, woo
            f=h/(M_PI/3.0)-(float)hp;

            p=v*(1.0-s);
            q=v*(1.0-f*s);
            t=v*(1.0-(1.0-f)*s);

            switch (hp){
                case 0: r=v; g=t; b=p; break;
                case 1: r=q; g=v; b=p; break;
                case 2: r=p; g=v; b=t; break;
                case 3: r=p; g=q; b=v; break;
                case 4: r=t; g=p; b=v; break;
                case 5: r=v; g=p; b=q; break;
            }

            //            fprintf(stderr,"h: %f r: %f g: %f b: %f\n",h,r,g,b);

  //          if (lattice[i][k][0].x == 0.0 && lattice[i][k][0].y == 0.0 && lattice[i][k][0].z == 0.0)
  //          { r=0.0; g=0.0; b=0.0; } // #FADE TO BLACK
            //zero length dipoles, i.e. absent ones - appear as black pixels

            fprintf(fo,"%c%c%c",(char)(254.0*r),(char)(254.0*g),(char)(254.0*b));
        }
    fclose(fo); //don't forget :^)
}

#define ZSCALE 5.0 // Scales Z-axis in Pymol xyz / CGO outputs

float DMAX=55.0; //sensible starting value...
float DMEAN=0.0;

void outputlattice_dumb_terminal()
{
    const char * species=".CZCT"; // Copper (I), Zinc (II), Tin (III)
    int x,y;
    float a;
    int z=0;
    float new_DMAX=0.0; //used to calibrate next colour scale, based on present maxima of data
    float variance=0.0; // sum of potential^2
    float mean=0.0;
    float potential;

    fprintf(stderr,"%*s%*s\n",X+3, "SPECIES", (2*X)+4,"POTENTIAL"); //padded labels


for (z=0;z<4;z++)
{
    fprintf(stderr,"Z=%d\n",z);

     for (y=0;y<Y;y++)
        for (x=0;x<X;x++)
        {
             potential=dipole_potential(x,y,z);
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
//            if (a<4.0)                                  // makes colour bold / normal depending on arrow orientation
                fprintf(stderr,";7"); // inverted colours
            char arrow=species[(int)a];
     //       if (lattice[x][y][z].z> sqrt(2)/2.0) arrow='o';
     //       if (lattice[x][y][z].z<-sqrt(2)/2.0) arrow='x';

   //         if (lattice[x][y][z].x==0.0 && lattice[x][y][z].y==0.0 && lattice[x][y][z].z==0.0) arrow='*'; 

            fprintf(stderr,"m%c %c[0m",arrow,27);  // prints arrow
            fprintf(stderr,"%c[37m%c[0m",27,27); //RESET

            //            fprintf(stderr,"%c ",arrows[(int)a]); // dumb - just black 'n'
            //            white
        }

        // OK - now potential plot :^)
        //        const char * density=".,:;o*O#"; //increasing potential density
        const char * density="012345689";
        fprintf(stderr,"    ");
        for (x=0;x<X;x++)
        {
            potential=dipole_potential(x,y,z);
            
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

            char arrow=species[(int)a];  // selectss arrow

            fprintf(stderr,"m%c%c%c[0m",density[(int)(8.0*fabs(potential-DMEAN)/DMAX)],arrow,27);
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
    int x,y,z;

    int dx=1,dy=1,dz=1; //which partner

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
            {
                printf("X:%d Y:%d Z:%d ",x,y,z);
                //int cutoff=3;
                for (int cutoff=1;cutoff<=7;cutoff++)
                {
                    double site_E=site_energy(x,y,z,lattice[x][y][z],cutoff);
                    printf("\n\t %d: %f",cutoff,site_E); 
                    
                    double dE=0.0; 
                    int species_a=lattice[x][y][z]; 
                    int species_b=lattice[x+dx][y+dy][z+dz];
                    
                    dE+=site_energy(x,y,z, species_a, cutoff);
                    dE-=site_energy(x,y,z, species_b, cutoff);
                    dE+=site_energy(x+dx,y+dy,z+dz, species_b, cutoff);
                    dE-=site_energy(x+dx,y+dy,z+dz, species_a, cutoff);
                    
                    printf(" %f",dE); 
                }
                printf("\n");
            }
   
}
