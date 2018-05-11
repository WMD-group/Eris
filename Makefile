SRCS=src/eris-kernel.c src/eris-analysis.c src/eris-config.c src/eris-lattice.c src/eris-main.c

eris: ${SRCS} 
	gcc -O4 -std=gnu11 -lm -lconfig -o eris src/eris-main.c

debug: ${SRCS}
	gcc -std=gnu11 -lm -lconfig -o eris -g src/eris-main.c

eris-openmp: ${SRCS} 
	gcc -O4 -std=gnu11 -lm -lconfig -fopenmp -o eris src/eris-main.c

eris-mac-openmp: ${SRCS} 
	/usr/local/bin/gcc-4.8 -std=gnu11 -O4 -lm -lconfig -fopenmp -lgomp -o eris src/eris-main.c

profile: ${SRCS} 
	gcc -std=gnu11 -lm -lconfig -o eris src/eris-main.c -pg

cx1-icc: 
	icc -std=c99 -Llibconfig-1.5/lib -Ilibconfig-1.5/lib \
	-O4 -o eris src/eris-main.c libconfig-1.5/lib/.libs/libconfig.a

parallel: eris 
	seq 0 50 1000 | parallel  ./eris {}

all: eris

clean:
	rm eris 
	
cleanupdata:
	rm czts* potential* Efield* variance* RDF*; rm -r equil* gulp_inputs

