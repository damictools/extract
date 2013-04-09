CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = extract.o 
HEADERS = globalConstants.h

ALL : extract.exe
	@echo "Listo!"

extract.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o extract.exe $(LIBS) $(GLIBS) $(CFLAGS)

extract.o : extract.cc $(HEADERS)
	$(CPP) -c extract.cc -o extract.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
