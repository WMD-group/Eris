/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#define X 30 // Malloc is for losers.
#define Y 30 
#define Z 10 

#define SPECIES 4

int DIM=3; //currently just whether the dipoles can point in Z-axis (still a 2D slab) 
int T; //global variable so accessible to analysis routines

struct dipole
{
    float x,y,z;
    float length; //length of dipole, to allow for solid state mixture (MA, FA, Ammonia, etc.)
}; 

int lattice[X][Y][Z];

double E_int[SPECIES][SPECIES]; // interaction energy between species

struct mixture
{
    float length;
    float prevalence;
} dipoles[10];
int dipolecount=0;

// SIMULATION PARAMETERS
// NB: These are defaults - most are now read from config file

double beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

struct dipole Efield; //now a vector, still k_B.T units per lattice unit
//double Efield=0.01; // units k_B.T per lattice unit
double Eangle=0.0;

double K=1.0; //elastic coupling constant for dipole moving within cage

double Dipole=1.0; //units of k_B.T for spacing = 1 lattice unit
double CageStrain=1.0; // as above

double dipole_fraction=0.9; //fraction of sites to be occupied by dipoles

int DipoleCutOff=1;

// These variables from the old main
int MCMegaSteps=400;
double MCMegaMultiplier=1.0;
int MCMinorSteps=0;
char const *LOGFILE = NULL; //for output filenames
//END OF SIMULATION PARAMETERS

// {{ Except for the ones hardcoded into the algorithm :^) }}

unsigned long ACCEPT=0; //counters for MC moves
unsigned long REJECT=0;

// Prototypes...
static float dot(struct dipole *a, struct dipole *b);
static void random_sphere_point(struct dipole *p);

// 3-Vector dot-product... hand coded, should probably validate against
// a proper linear albegra library
static float dot(struct dipole *a, struct dipole *b)
{
    int D;
    float sum=0.0;

    sum+=a->x*b->x;
    sum+=a->y*b->y;
    sum+=a->z*b->z;

    return(sum);
}

static void random_sphere_point(struct dipole *p)
{
    int i;
    // Marsaglia 1972 
    float x1,x2;
    do {
        x1=2.0*genrand_real1() - 1.0;
        x2=2.0*genrand_real1() - 1.0;
    } while (x1*x1 + x2*x2 > 1.0);

    if (DIM<3){
        // Circle picking, after Cook 1957
        // http://mathworld.wolfram.com/CirclePointPicking.html
        p->x = (x1*x1 - x2*x2)  / (x1*x1 + x2*x2);
        p->y =      2*x1*x2     / (x1*x1 + x2*x2);
        p->z = 0.0;
    }
    else
    {
        // Sphere picking
        p->x = 2*x1*sqrt(1-x1*x1-x2*x2);
        p->y = 2*x2*sqrt(1-x1*x1-x2*x2);
        p->z = 1.0 - 2.0* (x1*x1+x2*x2);
    }
}



void load_config()
{
    int i,j,k, x,y; //for loop iterators

    config_t cfg, *cf; //libconfig config structure
    const config_setting_t *setting;
    double tmp;
    int T;

// TODO: move this to init file

// Copper I
// Zinc II
// Tin (Sn) III
E_int[1][1]=-1.0;
E_int[1][2]=1.0;
E_int[1][3]=1.0;
E_int[2][1]=E_int[1][2];
E_int[2][2]=-1.0;
E_int[2][3]=1.0;
E_int[3][1]=E_int[1][3];
E_int[3][2]=E_int[2][3];
E_int[3][3]=-1.0;


    //Load and parse config file
    cf = &cfg;
    config_init(cf);

    if (!config_read_file(cf,"eris.cfg")) 
    {
        fprintf(stderr, "%s:%d - %s\n",
                config_error_file(cf),
                config_error_line(cf),
                config_error_text(cf));
        config_destroy(cf);
        exit(EXIT_FAILURE);
    }

    config_lookup_string(cf,"LOGFILE",&LOGFILE); //library does its own dynamic allocation

    config_lookup_int(cf,"T",&T);

    config_lookup_float(cf,"Efield.x",&tmp);  Efield.x=(float)tmp;
    config_lookup_float(cf,"Efield.y",&tmp);  Efield.y=(float)tmp;
    config_lookup_float(cf,"Efield.z",&tmp);  Efield.z=(float)tmp;

    fprintf(stderr,"Efield: x %f y %f z %f\n",Efield.x,Efield.y,Efield.z);

    //    config_lookup_float(cf,"Eangle",&Eangle);

    config_lookup_int(cf,"DipoleCutOff",&DipoleCutOff);

    config_lookup_int(cf,"MCMegaSteps",&MCMegaSteps);
    config_lookup_float(cf,"MCMegaMultiplier",&MCMegaMultiplier);

    MCMinorSteps=(int)((float)X*(float)Y*MCMegaMultiplier);

    fprintf(stderr,"Config loaded. \n");
}

