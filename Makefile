SRCS=eris-kernel.c eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c

eris: ${SRCS} 
	gcc -O4 -std=gnu11 -lm -lconfig -o eris eris-main.c

debug: ${SRCS}
	gcc -std=gnu11 -lm -lconfig -o eris -g eris-main.c

eris-openmp: ${SRCS} 
	gcc -O4 -std=gnu11 -lm -lconfig -fopenmp -o eris eris-main.c

eris-mac-openmp: ${SRCS} 
	/usr/local/bin/gcc-4.8 -std=gnu11 -O4 -lm -lconfig -fopenmp -lgomp -o eris eris-main.c

profile: ${SRCS} 
	gcc -std=gnu11 -lm -lconfig -o eris eris-main.c -pg

parallel: eris 
	seq 0 50 1000 | parallel  ./eris {}

all: eris

clean:
	rm eris 
	
cleandata:
	rm *.xyz potential*.dat variance.dat RDF_*.dat
