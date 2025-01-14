##	Makefile
#
#
### CLHEP
CLHEPDIR        = /usr/local
#
CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
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

OBJ =   ConfMan.o RawData.o \
	GetNumberFromKernelEntropyPool.o \
        PrimInfo.o HodoRawHit.o TrRawHit.o \
        s_BeamRawHit.o s_ScatRawHit.o \
        s_HodoRawHit.o s_TrRawHit.o \
        TrGeomMan.o TrGeomRecord.o \
        FieldMan.o SpecFieldMap.o FieldElements.o \
	Kinematics.o MathTools.o \
	Particle.o \
	HodoAnalyzer.o HodoHit.o HodoCluster.o\
	TrHit.o TrAnalyzer.o TrTrackSearch.o \
	TrHitCluster.o TrLocalTrack.o TrLTrackHit.o TrackHit.o \
        BeamTrack.o RungeKuttaUtilitiesBeam.o \
        PreInTrack.o RungeKuttaUtilitiesPreIn.o \
        PreOutTrack.o RungeKuttaUtilitiesPreOut.o \
        PreOut2Track.o RungeKuttaUtilitiesPreOut2.o \
        Scat1Track.o RungeKuttaUtilitiesScat1.o \
        Scat2Track.o RungeKuttaUtilitiesScat2.o \
        Scat2ATrack.o RungeKuttaUtilitiesScat2A.o \
        Scat2BTrack.o RungeKuttaUtilitiesScat2B.o \
        Scat3Track.o RungeKuttaUtilitiesScat3.o 

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysMonitor \
     analysLocalTracking \
     analysBeamTracking \
     analysScat1Tracking \
     analysScat2Tracking \
     analysScat3Tracking \
     analysSpecTracking \

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

analysMonitor: lib UserMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserMonitor.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysSMonitor: lib UserSMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSMonitor.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysLocalTracking: lib UserLocalTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserLocalTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysBeamTracking: lib UserBeamTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBeamTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysScat1Tracking: lib UserScat1Tracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserScat1Tracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysScat2Tracking: lib UserScat2Tracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserScat2Tracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysScat3Tracking: lib UserScat3Tracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserScat3Tracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysSpecTracking: lib UserSpecTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSpecTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libSpecAna.a ../bin/*
