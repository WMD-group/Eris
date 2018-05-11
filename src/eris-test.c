// eris-test.c
// Minimalist unit testing via the MinUnit framework:
// http://www.jera.com/techinfo/jtns/jtn002.html

/* file: minunit.h */
#define mu_assert(message, test) do { if (!(test)) return message; } while (0)
#define mu_run_test(test) do { char *message = test(); tests_run++; \
                               if (message) return message; } while (0)
extern int tests_run;
int tests_run=0;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libconfig.h>
#include "eris-random.c"
#include "eris-config.c" // TODO: bump everything into eris-kernel; except for bare main fn.

static char * test_dot() {
    struct dipole a,b;
    double dotted;

    a.x=1.0; a.y=0.0; a.z=0.0; a.length=1.0;
    b.x=1.0; b.y=0.0; b.z=0.0; b.length=1.0;

    dotted=dot(&a,&b);
    mu_assert("error, [1,0,0].[1,0,0] != 1", dotted == 1.0);
    
    return 0;
}

static char * all_tests() { 
	mu_run_test(test_dot);
    //mu_run_test(test_bar);
    return 0;
}
 
int main(int argc, char **argv) {
    char *result = all_tests();
    if (result != 0) {
        printf("%s\n", result);
    }
    else {
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);
 
    return result != 0;
}

