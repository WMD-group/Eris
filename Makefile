eris: eris-kernel.c eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c
	gcc -O4 -std=gnu11 -lm -lconfig -o eris eris-main.c

eris-openmp: eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c
	gcc -O4 -std=gnu11 -lm -lconfig -fopenmp -o eris eris-main.c

eris-mac-openmp: eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c
	/usr/local/bin/gcc-4.8 -std=gnu11 -O4 -lm -lconfig -fopenmp -lgomp -o eris eris-main.c

profile: eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c 
	gcc -std=gnu11 -lm -lconfig -o eris eris-main.c -pg

parallel: eris 
	seq 0 50 1000 | parallel  ./eris {}

all: eris

clean:
	rm eris 
	
cleanupdata:
	rm *.xyz potential*.dat variance.dat
