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

LIBOBJ = libAna.a
TESTEXLIB = -L. -lAna
LIBS = $(CLHEPLIBS) $(RTLIBS)
LIBSONL = $(CLHEPLIBS) $(RTLIBS) 

CC = gcc
C++ = g++

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o CMapMan.o GetNumberFromKernelEntropyPool.o \
        RawData.o Decoder.o HexDump.o \
        HodoRawHit.o HRTdcRawHit.o TrRawHit.o DCRawHit.o \
        TrGeomMan.o TrGeomRecord.o \
	DCHit.o DCAnalyzer.o DCLocalTrack.o DCLTrackHit.o \
	DCGeomMan.o DCGeomRecord.o DCTdcCalibMan.o DCDriftParamMan.o \
        DCLocalTrack.o DCLTrackHit.o TrackHit.o\
	DCTrackSearch.o DCPairHitCluster.o \
        MathTools.o

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysMonitor analysMonitor2 analysDCTracking

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

analysMonitor: lib UserMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserMonitor.o \
	Main.o $(TESTEXLIB) $(LIBS)

analysMonitor2: lib UserMonitor2.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserMonitor2.o \
	Main.o $(TESTEXLIB) $(LIBS)

analysDCTracking: lib UserDCTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserDCTracking.o \
	Main.o $(TESTEXLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libAna.a ../bin/*
