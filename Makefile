eris: eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c
	gcc -O4 -lm -lconfig -o eris eris-main.c

profile: eris-analysis.c   eris-config.c  eris-lattice.c  eris-main.c 
	gcc -lm -lconfig -o eris eris-main.c -pg

all: eris

clean:
	rm eris *.pnm *.jpg *.gif *.avi *.svg *.png
