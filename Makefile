CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags) 
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)  
GLIBS = 
GLIBS += 
OBJECTS = extract.o tinyxml2.o gConfig.o
HEADERS = globalConstants.h tinyxml2.h gConfig.h

ALL : extract.exe
	@echo "Listo!"

extract.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o extract.exe $(LIBS) $(GLIBS) $(CFLAGS)

extract.o : extract.cc $(HEADERS)
	$(CPP) -c extract.cc -o extract.o $(CFLAGS)

gConfig.o : gConfig.cc $(HEADERS)
	$(CPP) -c gConfig.cc -o gConfig.o $(CFLAGS)

tinyxml2.o : tinyxml2.cpp $(HEADERS)
	$(CPP) -c tinyxml2.cpp -o tinyxml2.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
