##	Makefile
#
#
### CLHEP
CLHEPDIR        = /cern/CLHEP
#
CLHEPINCDIR     = $(CLHEPDIR)/include
CLHEPLIBDIR     = $(CLHEPDIR)/lib
#
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIBS       = -L$(CLHEPLIBDIR) -lCLHEP
#
#
### RooT
RTINCLUDE	= `root-config --cflags`
RTLIBS		= `root-config --libs`
#
#
CFLAGS      = -g -O2 
#-Wall 
#DFLAGS      = -Df2cFortran -Dextname -DDEBUG 

#debug_flags += -DMemoryLeak

FLAGS = $(CFLAGS) $(DFLAGS) -I. $(CLHEPINCLUDE) $(RTINCLUDE)

LIBOBJ = libSpecAna.a
SPECLIB = -L. -lSpecAna
LIBS = $(CLHEPLIBS) $(RTLIBS)
LIBSONL = $(CLHEPLIBS) $(RTLIBS) 

CC = gcc
C++ = g++

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o \
        RawData.o HodoRawHit.o DCRawHit.o TrRawHit.o \
	GeomMan.o GeomRecord.o DCTdcCalibMan.o DCDriftParamMan.o \
	DCAnalyzer.o DCHit.o DCLocalTrack.o DCLTrackHit.o TrackHit.o \
	DCTrackSearch.o DCPairHitCluster.o \
	HodoAnalyzer.o HodoParamMan.o HodoPHCMan.o \
	Hodo2Hit.o Hodo1Hit.o HodoCluster.o \
	MathTools.o Kinematics.o

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysMonitor analysHodoscope \
     analysKLDCTracking analysBDCTracking \
     analysReacAna01

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

analysMonitor: lib UserMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserMonitor.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysHodoscope: lib UserHodoscope.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserHodoscope.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysKLDCTracking: lib UserKLDCTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserKLDCTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysBDCTracking: lib UserBDCTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBDCTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysReacAna01: lib UserReacAna01.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserReacAna01.o 	\
	Main.o $(SPECLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libSpecAna.a ../bin/*
