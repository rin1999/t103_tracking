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
### HD DAQ Unpacker ####
unpacker_install      := $(HOME)/K1.8Analyzer/unpacker
unpacker_include_dir  := $(unpacker_install)/include
unpacker_cflags       := -I$(unpacker_include_dir)
unpacker_lib_dir      := $(unpacker_install)/lib
unpacker_libs         := \
 -Wl,-rpath,$(unpacker_lib_dir) -L$(unpacker_lib_dir) -lHDDAQUnpacker
#
#
CFLAGS      = -g -O2 
#-Wall 
#DFLAGS      = -Df2cFortran -Dextname -DDEBUG 

#debug_flags += -DMemoryLeak

FLAGS = $(CFLAGS) $(DFLAGS) -I. $(CLHEPINCLUDE) $(RTINCLUDE)\
	$(unpacker_cflags) 
#$(debug_flags)

LIBOBJ = libSksAna.a
SKSLIB = -L. -lSksAna
LIBS = $(CLHEPLIBS) $(RTLIBS) $(unpacker_libs)
LIBSONL = $(CLHEPLIBS) $(RTLIBS) $(unpacker_libs)

C++ = g++

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o RawData.o \
	HodoRawHit.o HodoAnalyzer.o HodoParamMan.o HodoPHCMan.o \
	Hodo2Hit.o Hodo1Hit.o\
	BH2Hit.o \
	HodoCluster.o BH2Cluster.o \
	DCRawHit.o DCHit.o DCAnalyzer.o DCLocalTrack.o DCLTrackHit.o\
	TrackHit.o MWPCCluster.o TrackMaker.o BH1Filter.o \
	K18Track.o K18TrackFCN.o K18TransMatrix.o K18Particle.o Minuit.o\
	SksTrack.o SksParticle.o RungeKuttaUtilities.o MathTools.o\
	FieldMan.o SksFieldMap.o FieldElements.o \
	DCGeomMan.o DCGeomRecord.o DCTdcCalibMan.o DCDriftParamMan.o \
	DCTrackSearch.o DCPairHitCluster.o\
	Kinematics.o DCTrackUtilities.o EvDisp.o \
        DebugTimer.o DebugCounter.o
# 	ScalerAna.o\
# 	DSio.o CFile.o EvDisp.o   		\
# 	Sdc1HitCluster.o  D           \
# 	SimuData.o 	\
# 	

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysMonitor analysEvDisp \
     analysHodoMonitor analysHodoscope \
     analysDCMonitor analysDCCheck analysDCHitCheck \
     analysSdcInTracking analysSdcOutTracking \
     analysSksTracking analysSksHodoscope\
     analysSksHodoscopeTheta analysSksHodoscopePiP \
     analysSksHodoscopeSigma analysSksHodoscopeSigmaCut \
     analysBcSdcBeamThrough \
     analysBcInTracking analysBcOutTracking \
     analysK18Tracking \
     analysBeamThrough \
     analysPiPAna analysPiPiAna analysPiKAna \
     analysPiKAnaSigma \
     analysPiKAnaSigmaPlus analysPiKAnaSigmaMinus \
     analysPiKAnaTheta \
     analysPiPiAnaTheta \
     analysHodoscopeSP0 	

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

analysMonitor: lib UserMonitor.o Main.o $(LIBOBJ)
	$(C++) -o $(BINDIR)/$@ $(FLAGS) UserMonitor.o Main.o \
	$(SKSLIB) $(LIBS) 

analysEvDisp: lib UserEvDisp.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserEvDisp.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysHodoMonitor: lib UserHodoMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserHodoMonitor.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysHodoscope: lib UserHodoscope.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserHodoscope.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysHodoscopeSP0: lib UserHodoscopeSP0.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserHodoscopeSP0.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysDCMonitor: lib UserDCMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserDCMonitor.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysDCCheck: lib UserDCCheck.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserDCCheck.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysDCHitCheck: lib UserDCHitCheck.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserDCHitCheck.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSdcInTracking: lib UserSdcInTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSdcInTracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSdcOutTracking: lib UserSdcOutTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSdcOutTracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksTracking: lib UserSksTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksTracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksHodoscope: lib UserSksHodoscope.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksHodoscope.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksHodoscopeTheta: lib UserSksHodoscopeTheta.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksHodoscopeTheta.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksHodoscopeSigma: lib UserSksHodoscopeSigma.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksHodoscopeSigma.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksHodoscopeSigmaCut: lib UserSksHodoscopeSigmaCut.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksHodoscopeSigmaCut.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysSksHodoscopePiP: lib UserSksHodoscopePiP.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSksHodoscopePiP.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysBcInTracking: lib UserBcInTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBcInTracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysBcOutTracking: lib UserBcOutTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBcOutTracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysK18Tracking: lib UserK18Tracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserK18Tracking.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysBcSdcBeamThrough: lib UserBcSdcBeamThrough.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBcSdcBeamThrough.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysBeamThrough: lib UserBeamThrough.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBeamThrough.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiPAna: lib UserPiPAna.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiPAna.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiPiAna: lib UserPiPiAna.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiPiAna.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiKAna: lib UserPiKAna.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiKAna.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiKAnaSigma: lib UserPiKAnaSigma.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiKAnaSigma.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiKAnaSigmaPlus: lib UserPiKAnaSigmaPlus.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiKAnaSigmaPlus.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiKAnaSigmaMinus: lib UserPiKAnaSigmaMinus.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiKAnaSigmaMinus.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiKAnaTheta: lib UserPiKAnaTheta.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiKAnaTheta.o 	\
	Main.o $(SKSLIB) $(LIBS)

analysPiPiAnaTheta: lib UserPiPiAnaTheta.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPiPiAnaTheta.o 	\
	Main.o $(SKSLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libSksAna.a
