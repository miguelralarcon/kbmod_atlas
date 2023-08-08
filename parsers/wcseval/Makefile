CC= cc
CFLAGS= -O -g -Wall

pix2sky: pix2sky.o wcseval.o
	$(CC) $(CFLAGS) pix2sky.o wcseval.o -lm -o pix2sky

install:
	install -m 0775 pix2sky $(ATLAS_HOME)/bin/
	ln -fs $(ATLAS_HOME)/bin/pix2sky $(ATLAS_HOME)/bin/sky2pix
	install -m 0664 wcseval.o $(ATLAS_HOME)/lib/
	install -m 0664 wcseval.h $(ATLAS_HOME)/include/

all: pix2sky wcseval.o

clean:
	rm -f *.o pix2sky
