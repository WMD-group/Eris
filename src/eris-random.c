/* eris-random.c
 * Wrapped for compiled random number generator + convenience functions.
 */

#include "mt19937ar-cok.c" //Code _included_ to allow more global optimisation
static int rand_int(int SPAN) // TODO: profile this to make sure it runs at an OK speed.
{
        return((int)( (unsigned long) genrand_int32() % (unsigned long)SPAN));
}

