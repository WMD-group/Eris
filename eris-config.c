/* Starry Night - a Monte Carlo code to simulate ferroelectric domain formation
 * and behaviour in hybrid perovskite solar cells.
 *
 * By Jarvist Moore Frost
 * University of Bath
 *
 * File begun 16th January 2014
 */

#include <stdbool.h>

// ----------------------------------------------------------------------------------------------------
// SIMULATION PARAMETERS
// NB: These are defaults - most are now read from config file

double beta=1.0;  // beta=1/T  T=temperature of the lattice, in units of k_B

int ElectrostaticCutOff=1;
int freezeSn = true; // should we freeze Sn ?

int MCMegaSteps=400;
int MCEqmSteps=0;
int TMAX=500;
int TMIN=0;
int TSTEP=100; 

double MCMoves=1.0;
unsigned long long int MCMinorSteps=0;
char const *LOGFILE = NULL; //for output filenames

int DEBUG=false;
int DisplayDumbTerminal=true;
int DumbTerminalLayers=4;
int CalculateRadialOrderParameter=false;
int CalculatePotential=false;
int OrderedInitialLattice=false;
int ReinitialiseLattice=false;
int EquilibrationChecks=false;

int SuzySupercell=false;

int SaveXYZ=false;
int SaveGULP=false;
//END OF SIMULATION PARAMETERS
// ----------------------------------------------------------------------------------------------------------------------------------


// User defined lattice dimensions for both original lattice initialisation method and for SuzySupercell method ---------------------- 
// ----------------------------------------------------------------------------------------------------------------------------------

void set_lattice_dimensions(int * X, Y, Z);

void set_lattice_dimensions(int * X, Y, Z)
{

// Code for reading eris.cfg (just for setting lattice dimensions), file read fully in load_config function later-------------------
// ---------------------------------------------------------------------------------------------------------------------------------
  config_t cfg, *cf; //libconfig config structure
  const config_setting_t *setting;
    
  //Load and parse config file
  cf = &cfg;
  config_init(cf);

  fprintf(stderr,"Reading config file...\n");
  if (!config_read_file(cf,"eris.cfg"))
  {
      fprintf(stderr, "%s:%d - %s\n",
              config_error_file(cf),
              config_error_line(cf),
              config_error_text(cf));
      config_destroy(cf);
      exit(EXIT_FAILURE);
  }
// -------------------------------------------------------------------------------------

  // Checking eris.cfg for which lattice initialisation method the user has selected (SuzySupercell or original method if false)
  config_lookup_bool(cf,"SuzySupercell",&SuzySupercell);

  if (SuzySupercell)
  {
    // User defined system dimensions to create a supercell of a 2x2x4 unit cell using SuzySupercell lattice initialisation method
    #define X_super 2  // set supercell dimension here
    #define Y_super 2  // and here
    #define Z_super 2  // and here
    int X = X_super*2;
    int Y = Y_super*2;
    int Z = Z_super*4;
    int lattice[X][Y][Z];
  }  
  else
  {
    // User defined dimensions if using original lattice initialisation method
    #define X 10 // Malloc is for losers.
    #define Y 10 // X must be divisible by 4, Y divisible by 2, to generate stoichometric CZTS 
    #define Z 20 
    int lattice[X][Y][Z];
  }
}



#define POTENTIAL_CUTOFF 4 // cutoff for calculation of electrostatic potential
// ill defined if this is > than half any of the above

#define SPECIES 4

int DIM=3; //currently just whether the dipoles can point in Z-axis (still a 2D slab) 
int T; //global variable so accessible to analysis routines

// These defines set enum values for the species type; so that we don't need
// magic numbers.
// BUT - it makes Eris specific to CZTS material.
#define Cu 1
#define Zn 2
#define Sn 3
// And similarly for the reverse lookup (single character), i.e. species[2]='Z'
const char * specieslookup=".CZTc"; // Copper = 1, Zinc = 2, Tin = 3

struct dipole
{
    float x,y,z;
    float length; //length of dipole, to allow for solid state mixture (MA, FA, Ammonia, etc.)
}; 

double E_int[SPECIES][SPECIES]; // interaction energy between species
double FormalCharge[SPECIES];


struct mixture
{
    float length;
    float prevalence;
} dipoles[10];
int dipolecount=0;


// {{ Except for the ones hardcoded into the algorithm :^) }}

unsigned long long int ACCEPT=0; //counters for MC moves
unsigned long long int REJECT=0;

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


// Load eris.cfg; overwriting global variable defaults above.
void load_config()
{
    int i,j,k, x,y; //for loop iterators

    config_t cfg, *cf; //libconfig config structure
    const config_setting_t *setting;
    int E_ints;
    double electrostatic;
    double tmp;

    //Load and parse config file
    cf = &cfg;
    config_init(cf);

    fprintf(stderr,"Reading config file...\n");
    if (!config_read_file(cf,"eris.cfg")) 
    {
        fprintf(stderr, "%s:%d - %s\n",
                config_error_file(cf),
                config_error_line(cf),
                config_error_text(cf));
        config_destroy(cf);
        exit(EXIT_FAILURE);
    }
    fprintf(stderr,"Read OK! Variables have been set as follows:\n\n");

    config_lookup_string(cf,"LOGFILE",&LOGFILE); //config library does its own dynamic allocation
    config_lookup_int(cf,"T",&T); 
    config_lookup_float(cf,"electrostatic",&electrostatic); //Multiplier for E_ints
    fprintf(stderr,"LOGFILE: %s\nT: %d\nelectrostatic: %f\n",LOGFILE,T,electrostatic);

    setting = config_lookup(cf, "E_int");
    E_ints   = config_setting_length(setting);
    fprintf(stderr,"I've found: %d values in the E_int list...\n",E_ints);
//   Nb: factor of 38.911 kBT at 300 K in eV; to convert eV values into internal ones...
    for (i=0;i<E_ints;i++)
        E_int[i/3][i%3]=config_setting_get_float_elem(setting,i)*electrostatic*38.911; //I know, I know - I'm sorry.
    //    config_lookup_float(cf,"Eangle",&Eangle);
    fprintf(stderr,"My interaction matrix (should be symmetric) look like:\n");
    for (i=0;i<E_ints/3;i++)
        fprintf(stderr,"    %f %f %f\n",E_int[i][0],E_int[i][1],E_int[i][2]);

    int FormalCharges;
    setting  = config_lookup(cf,"FormalCharges");
    FormalCharges = config_setting_length(setting);
    fprintf(stderr,"I've found: %d values in the FormalCharges list...\n",FormalCharges);
    for (i=0;i<FormalCharges;i++)
        FormalCharge[i]=config_setting_get_float_elem(setting,i);
    fprintf(stderr,"Formal charges look like: ");
    for (i=0;i<FormalCharges;i++)
        fprintf(stderr,"%f\t",FormalCharge[i]);
    fprintf(stderr,"\n");

    config_lookup_int(cf,"ElectrostaticCutOff",&ElectrostaticCutOff);
    config_lookup_int(cf,"MCMegaSteps",&MCMegaSteps);
    config_lookup_int(cf,"MCEqmSteps",&MCEqmSteps);
    config_lookup_float(cf,"MCMoves",&MCMoves);
    fprintf(stderr,"ElectrostaticCutOff: %d\n\nMCMegaSteps: %d\nMCEqmSteps: %d\nMCMoves: %f\n",ElectrostaticCutOff,MCMegaSteps,MCEqmSteps,MCMoves);

    MCMinorSteps=(unsigned long long int)((float)X*(float)Y*(float)Z*MCMoves);
    fprintf(stderr," ==> Results in MCMinorSteps: %llu\n",MCMinorSteps);

    config_lookup_int(cf,"TMIN",&TMIN);
    config_lookup_int(cf,"TMAX",&TMAX);
    config_lookup_int(cf,"TSTEP",&TSTEP);
    fprintf(stderr,"\nTMIN: %d\nTMAX: %d\nTSTEP: %d\n",TMIN,TMAX,TSTEP);
    if (TMIN<0 || TMAX<0 || TSTEP<1)
        fprintf(stderr,"SOMETHING VERY ODD ABOUT THE TEMPERATURES I READ FROM THE CONFIG FILE. I HOPE YOU KNOW WHAT YOU ARE DOING!\n");

// Flags for output routines to run
    config_lookup_bool(cf,"DEBUG",&DEBUG);
    config_lookup_bool(cf,"DisplayDumbTerminal",&DisplayDumbTerminal);
    config_lookup_int(cf,"DumbTerminalLayers",&DumbTerminalLayers); 
    config_lookup_bool(cf,"CalculateRadialOrderParameter",&CalculateRadialOrderParameter);
    config_lookup_bool(cf,"CalculatePotential",&CalculatePotential);
    config_lookup_bool(cf,"OrderedInitialLattice",&OrderedInitialLattice);
    config_lookup_bool(cf,"ReinitialiseLattice",&ReinitialiseLattice);
    
    config_lookup_bool(cf, "EquilibrationChecks",&EquilibrationChecks);

    config_lookup_bool(cf,"SuzySupercell",&SuzySupercell);

    config_lookup_bool(cf,"SaveXYZ",&SaveXYZ);
    config_lookup_bool(cf,"SaveGULP",&SaveGULP);

    config_lookup_bool(cf,"freezeSn",&freezeSn);

    fprintf(stderr,"Config loaded. \n\n");
}

