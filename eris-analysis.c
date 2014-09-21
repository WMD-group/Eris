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
    int MAX=10;
    double pot=0.0;
    float d;
    struct dipole r;

    for (dx=-MAX;dx<MAX;dx++)
        for (dy=-MAX;dy<MAX;dy++)
#if(Z>1) //i.e. 3D in Z
            for (dz=-MAX;dz<MAX;dz++)
#endif
            {
                if (dx==0 && dy==0 && dz==0)
                    continue; //no infinities / self interactions please!

                r.x=(float)(dx); r.y=(float)(dy); r.z=(float)(dz);

                d=sqrt((float) r.x*r.x + r.y*r.y + r.z*r.z); //that old chestnut

                if (d>(float)MAX) continue; // Cutoff in d

                // pot(r) = 1/4PiEpsilon * p.r / r^3
                // Electric dipole potential
// FIXME
                //                pot+=dot(& lattice[(X+x+dx)%X][(Y+y+dy)%Y][(Z+z+dz)%Z] ,& r)/(d*d*d);
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
    double pot;
    FILE *fo;
    fo=fopen(filename,"w");

    for (x=0;x<X;x++)
        for (y=0;y<Y;y++)
            for (z=0;z<Z;z++)
                fprintf(fo,"%d %d %f\n",x,y,dipole_potential(x,y,z));
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

void outputlattice_dumb_terminal()
{
    const char * species=".CZT"; // Copper (I), Zinc (II), Tin (III)
    int x,y;
    float a;
    int z=0;
    float new_DMAX=0.0; //used to calibrate next colour scale, based on present maxima of data
    float variance=0.0; // sum of potential^2
    float mean=0.0;

    fprintf(stderr,"%*s%*s\n",X+3, "DIPOLES", (2*X)+4,"POTENTIAL"); //padded labels

    for (y=0;y<Y;y++)
    {
        for (x=0;x<X;x++)
        {
            a=lattice[x][y][0];

            fprintf (stderr,"%c[%d",27,31+((int)a)%8 ); // Sets colour of output routine
            if (a<4.0)                                  // makes colour bold / normal depending on arrow orientation
                fprintf(stderr,";7");
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
        float potential;
        //        const char * density=".,:;o*O#"; //increasing potential density
        const char * density="012345689";
        fprintf(stderr,"    ");
        for (x=0;x<X;x++)
        {
            potential=dipole_potential(x,y,z);
            variance+=potential*potential;
            mean+=potential;

            if (fabs(potential)>new_DMAX)
                new_DMAX=fabs(potential); // used to calibrate scale - technically this changes
            //printf("%f\t",potential); //debug routine to get scale

            //fprintf(stderr,"%c[%d",27,31+((int)(8.0*fabs(potential)/DMAX))%8); //8 colours
            //fprintf(stderr,"%c[48;5;%d",27,17+(int)(214.0*fabs(potential)/DMAX)); // Xterm 256 color map - (16..231)
            fprintf(stderr,"%c[48;5;%d",27,232+12+(int)(12.0*potential/DMAX)); // Xterm 256 color map - shades of grey (232..255)
            // https://code.google.com/p/conemu-maximus5/wiki/AnsiEscapeCodes#xterm_256_color_processing_requirements

            //if (potential<0.0) // if negative
            //    fprintf(stderr,";7"); // bold

            a=lattice[x][y][0];

            char arrow=species[(int)a];  // selectss arrow
            //if (lattice[x][y][z].z> sqrt(2)/2.0) arrow='o'; // override for 'up' (towards you - physics arrow style 'o')
            //if (lattice[x][y][z].z<-sqrt(2)/2.0) arrow='x'; // and 'down' (away from you, physics arrow style 'x')


            fprintf(stderr,"m%c%c%c[0m",density[(int)(8.0*fabs(potential)/DMAX)],arrow,27);
        }

        fprintf(stderr,"\n");
    }
    mean=mean/(X*Y);
    variance=variance/(X*Y); 
    fprintf(stderr,"T: %d DMAX: %f new_DMAX: %f variance: %f mean: %f\n",T,DMAX,new_DMAX,variance,mean);
    fprintf(stdout,"T: %d DMAX: %f new_DMAX: %f variance: %f mean: %f\n",T,DMAX,new_DMAX,variance,mean);
    DMAX=(new_DMAX+DMAX)/2.0; // mean of old and new (sampled, noisy) value
    DMAX=new_DMAX; // infinite fast following - but leads to fluctuations at steady state
    if (DMAX==0.0) DMAX=1.0; //avoid divide by zero for all-zero pot
}
