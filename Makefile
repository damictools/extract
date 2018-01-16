CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS =
GLIBS +=
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++11
else
  CFLAGS += -std=c++0x
endif
OBJECTS = extract.o tinyxml2.o gConfig.o globalConstants.o 
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

globalConstants.o : globalConstants.cc $(HEADERS)
	$(CPP) -c globalConstants.cc -o globalConstants.o $(CFLAGS)
	
clean:
	rm -f *~ *.o *.exe
