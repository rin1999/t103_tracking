/*
  UserSMonitor.cc

  2016/2 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "PrimInfo.hh"
#include "HodoRawHit.hh"
#include "TrRawHit.hh"
#include "s_BeamRawHit.hh"
#include "s_ScatRawHit.hh"
#include "s_HodoRawHit.hh"
#include "s_TrRawHit.hh"

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventSMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventSMonitor();
  ~EventSMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventSMonitor::EventSMonitor()
  : VEvent(),
    rawData(0)
{
}

EventSMonitor::~EventSMonitor()
{
  if (rawData) delete rawData;
}

struct Event{
  //Primary
  int    prinhits;
  double priposx, priposy, priposz;
  double pbeam, ubeam, vbeam, abmom;
  double prim1, prip1;
  double pritheta1, priphi1, prithetacm1, priphicm1;
  double prim2, prip2;
  double pritheta2, priphi2, prithetacm2, priphicm2;

  //Generated Beam
  int    gbnhits;
  std::vector<int>    gbid, gbtype;
  std::vector<double> gbp, gbpx, gbpy, gbpz;
  std::vector<int>    gbpid;
  std::vector<double> gbvx, gbvy;

  //Generated Scattered
  int    gsnhits;
  std::vector<int>    gsid, gstype;
  std::vector<double> gsp, gspx, gspy, gspz;
  std::vector<int>    gspid;
  std::vector<double> gsvx, gsvy;

  //BFT
  int    bftnhits;
  std::vector<int> bftlayer;
  std::vector<double> bftposx, bftposy;
  std::vector<double> bftdl;

  //SFT
  int    sftnhits;
  std::vector<int> sftlayer;
  std::vector<double> sftposx, sftposy;
  std::vector<double> sftdl;

  //AFT
  int    aftnhits;
  std::vector<int> aftlayer;
  std::vector<double> aftposx, aftposy;
  std::vector<double> aftdl;

  //IT1
  int    it1nhits;
  std::vector<int> it1layer;
  std::vector<double> it1posx, it1posy;
  std::vector<double> it1dl;

  //IT2
  int    it2nhits;
  std::vector<int> it2layer;
  std::vector<double> it2posx, it2posy;
  std::vector<double> it2dl;

  //ST1
  int    st1nhits;
  std::vector<int> st1layer;
  std::vector<double> st1posx, st1posy;
  std::vector<double> st1dl;

  //ST2
  int    st2nhits;
  std::vector<int> st2layer;
  std::vector<double> st2posx, st2posy;
  std::vector<double> st2dl;

  //T0
  int    t0nhits;
  std::vector<int>    t0layer, t0seg;
  std::vector<double> t0time, t0edep;
  std::vector<double> t0path, t0p;
  std::vector<double> t0posx, t0posy;
  std::vector<int>    t0pid;
  std::vector<double> t0beta;

  //TOF
  int    tofnhits;
  std::vector<int>    toflayer, tofseg;
  std::vector<double> toftime, tofedep;
  std::vector<double> tofpath, tofp;
  std::vector<double> tofposx, tofposy;
  std::vector<int>    tofpid;
  std::vector<double> tofbeta;

  //ITOF
  int    itofnhits;
  std::vector<int>    itoflayer, itofseg;
  std::vector<double> itoftime, itofedep;
  std::vector<double> itofpath, itofp;
  std::vector<double> itofposx, itofposy;
  std::vector<int>    itofpid;
  std::vector<double> itofbeta;

  //PAD
  int    padnhits;
  std::vector<int>    padlayer, padseg;
  std::vector<double> padtime, padedep;
  std::vector<double> padpath, padp;
  std::vector<double> padposx, padposy;
  std::vector<int>    padpid;
  std::vector<double> padbeta;

  //RICH
  int    richnhits;
  std::vector<int>    richlayer, richseg;
  std::vector<double> richtime, richedep;
  std::vector<double> richpath, richp;
  std::vector<double> richposx, richposy;
  std::vector<int>    richpid;
  std::vector<double> richbeta;

  //PID1
  int    pid1nhits;
  std::vector<int>    pid1layer, pid1seg;
  std::vector<double> pid1time, pid1edep;
  std::vector<double> pid1path, pid1p;
  std::vector<double> pid1posx, pid1posy;
  std::vector<int>    pid1pid;
  std::vector<double> pid1beta;

  //PID2
  int    pid2nhits;
  std::vector<int>    pid2layer, pid2seg;
  std::vector<double> pid2time, pid2edep;
  std::vector<double> pid2path, pid2p;
  std::vector<double> pid2posx, pid2posy;
  std::vector<int>    pid2pid;
  std::vector<double> pid2beta;

  //MF
  int    mfnhits;
  std::vector<int>    mflayer, mfseg;
  std::vector<double> mftime, mfedep;
  std::vector<double> mfpath, mfp;
  std::vector<double> mfposx, mfposy;
  std::vector<int>    mfpid;
  std::vector<double> mfbeta;

  //VD
  int    vdnhits;
  std::vector<int>    vdlayer, vdseg;
  std::vector<double> vdtime, vdedep;
  std::vector<double> vdpath, vdp;
  std::vector<double> vdposx, vdposy;
  std::vector<int>    vdpid;
  std::vector<double> vdbeta;
};
static Event event;

bool EventSMonitor::ProcessingBegin()
{
 return true;
}

bool EventSMonitor::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits(In) ) return false;
  //std::cout << "***" << std::endl;

  //**************************************************************************
  //******************RawData

  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //PrimaryInfo
  {
    const PrimInfoContainer &cont=rawData->GetPrimHC();  
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      PrimInfo *hit=cont[i];
      event.priposx = hit->GetVertX();
      event.priposy = hit->GetVertY();
      event.priposz = hit->GetVertZ();
      
      event.pbeam = hit->GetBeamMom();
      event.ubeam = hit->GetBeamU();
      event.vbeam = hit->GetBeamV();
      event.abmom = hit->GetAnaBeamMom();
      
      event.prim1 = hit->GetMass1();
      event.prip1 = hit->GetMom1();
      event.pritheta1 = hit->GetTheta1();
      event.priphi1 = hit->GetPhi1();
      event.prithetacm1 = hit->GetThetaCM1();
      event.priphicm1 = hit->GetPhiCM1();
      
      event.prim2 = hit->GetMass2();
      event.prip2 = hit->GetMom2();
      event.pritheta2 = hit->GetTheta2();
      event.priphi2 = hit->GetPhi2();
      event.prithetacm2 = hit->GetThetaCM2();
      event.priphicm2 = hit->GetPhiCM2();
    }
  }  

  //Beam
  {
    const s_BeamRHitContainer &cont =rawData->GetsBeamRHC();
    int nhB = cont.size();
    event.gbnhits = nhB;
    for(int i=0; i<nhB; i++){
      s_BeamRawHit *hit=cont[i];
      event.gbid.push_back(hit->TrackId());
      event.gbtype.push_back(hit->TrackType());
      event.gbp.push_back(hit->GetMom());
      event.gbpx.push_back(hit->GetMomX());
      event.gbpy.push_back(hit->GetMomY());
      event.gbpz.push_back(hit->GetMomZ());
      event.gbpid.push_back(hit->GetPid());
      event.gbvx.push_back(hit->GetVertX());
      event.gbvy.push_back(hit->GetVertY());

      //T0
      const s_HodoRHitContainer &contT0 =hit->GetsT0RHC();
      int nhT0 = contT0.size();
      event.t0nhits = nhT0;
      for(int j=0; j<nhT0; j++){
      	s_HodoRawHit *hitT0=contT0[i];
      	event.t0layer.push_back(hitT0->LayerId());
      	event.t0seg.push_back(hitT0->SegmentId());

      	int nt = hitT0->GetSize();
      	for( int k=0; k<nt; k++ ){
      	  event.t0time.push_back(hitT0->GetTime(k));
      	  event.t0edep.push_back(hitT0->GetEdep(k));
      	  event.t0path.push_back(hitT0->GetPath(k));
      	  event.t0p.push_back(hitT0->GetMom(k));
      	  event.t0posx.push_back(hitT0->GetPosX(k));
      	  event.t0posy.push_back(hitT0->GetPosY(k));
      	  event.t0pid.push_back(hitT0->GetPid(k));
      	  event.t0beta.push_back(hitT0->GetBeta(k));
      	}
      }

      //BFT
      for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
      	const s_TrRHitContainer &contBFT = hit->GetsBFTRHC(layer);
      	int nhBFT=cont.size();
      	event.bftnhits = nhBFT;
      	for( int j=0; j<nhBFT; ++j ){
      	  s_TrRawHit *hitBFT=contBFT[j];
      	  event.bftlayer.push_back(hitBFT->LayerId());

      	  int nt = hitBFT->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.bftposx.push_back(hitBFT->GetPosX(k));
      	    event.bftposy.push_back(hitBFT->GetPosY(k));
      	    event.bftdl.push_back(hitBFT->GetDL(k));
      	  }
      	}
      }
    }
  }

  //Scattered
  {
    const s_ScatRHitContainer &cont =rawData->GetsScatRHC();
    int nhS = cont.size();
    event.gsnhits = nhS;
    for(int i=0; i<nhS; i++){
      s_ScatRawHit *hit=cont[i];

      event.gsid.push_back(hit->TrackId());
      event.gstype.push_back(hit->TrackType());
      event.gsp.push_back(hit->GetMom());
      event.gspx.push_back(hit->GetMomX());
      event.gspy.push_back(hit->GetMomY());
      event.gspz.push_back(hit->GetMomZ());
      event.gspid.push_back(hit->GetPid());
      event.gsvx.push_back(hit->GetVertX());
      event.gsvy.push_back(hit->GetVertY());

      //std::cout<< "*****" << std::endl;
      //std::cout<< "id= " << hit->TrackId() << std::endl;

      //TOF
      const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
      int nhTOF = contTOF.size();
      event.tofnhits = nhTOF;
      for(int j=0; j<nhTOF; j++){
       	s_HodoRawHit *hitTOF=contTOF[j];
	event.toflayer.push_back(hitTOF->LayerId());
      	event.tofseg.push_back(hitTOF->SegmentId());

	//std::cout<< "toflayer= " << hitTOF->LayerId() << std::endl;

      	int nt = hitTOF->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.toftime.push_back(hitTOF->GetTime(k));
      	  event.tofedep.push_back(hitTOF->GetEdep(k));
      	  event.tofpath.push_back(hitTOF->GetPath(k));
      	  event.tofp.push_back(hitTOF->GetMom(k));
      	  event.tofposx.push_back(hitTOF->GetPosX(k));
      	  event.tofposy.push_back(hitTOF->GetPosY(k));
      	  event.tofpid.push_back(hitTOF->GetPid(k));
      	  event.tofbeta.push_back(hitTOF->GetBeta(k));
	}
      }

      //ITOF
      const s_HodoRHitContainer &contITOF =hit->GetsITOFRHC();
      int nhITOF = contITOF.size();
      event.itofnhits = nhITOF;
      for(int j=0; j<nhITOF; j++){
       	s_HodoRawHit *hitITOF=contITOF[j];
	event.itoflayer.push_back(hitITOF->LayerId());
      	event.itofseg.push_back(hitITOF->SegmentId());

      	int nt = hitITOF->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.itoftime.push_back(hitITOF->GetTime(k));
      	  event.itofedep.push_back(hitITOF->GetEdep(k));
      	  event.itofpath.push_back(hitITOF->GetPath(k));
      	  event.itofp.push_back(hitITOF->GetMom(k));
      	  event.itofposx.push_back(hitITOF->GetPosX(k));
      	  event.itofposy.push_back(hitITOF->GetPosY(k));
      	  event.itofpid.push_back(hitITOF->GetPid(k));
      	  event.itofbeta.push_back(hitITOF->GetBeta(k));
	}
      }

      //PAD
      const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
      int nhPAD = contPAD.size();
      event.padnhits = nhPAD;
      for(int j=0; j<nhPAD; j++){
       	s_HodoRawHit *hitPAD=contPAD[j];
	event.padlayer.push_back(hitPAD->LayerId());
      	event.padseg.push_back(hitPAD->SegmentId());

      	int nt = hitPAD->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.padtime.push_back(hitPAD->GetTime(k));
      	  event.padedep.push_back(hitPAD->GetEdep(k));
      	  event.padpath.push_back(hitPAD->GetPath(k));
      	  event.padp.push_back(hitPAD->GetMom(k));
      	  event.padposx.push_back(hitPAD->GetPosX(k));
      	  event.padposy.push_back(hitPAD->GetPosY(k));
      	  event.padpid.push_back(hitPAD->GetPid(k));
      	  event.padbeta.push_back(hitPAD->GetBeta(k));
	}
      }

      //RICH
      const s_HodoRHitContainer &contRICH =hit->GetsRICHRHC();
      int nhRICH = contRICH.size();
      event.richnhits = nhRICH;
      for(int j=0; j<nhRICH; j++){
       	s_HodoRawHit *hitRICH=contRICH[j];
	event.richlayer.push_back(hitRICH->LayerId());
      	event.richseg.push_back(hitRICH->SegmentId());

      	int nt = hitRICH->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.richtime.push_back(hitRICH->GetTime(k));
      	  event.richedep.push_back(hitRICH->GetEdep(k));
      	  event.richpath.push_back(hitRICH->GetPath(k));
      	  event.richp.push_back(hitRICH->GetMom(k));
      	  event.richposx.push_back(hitRICH->GetPosX(k));
      	  event.richposy.push_back(hitRICH->GetPosY(k));
      	  event.richpid.push_back(hitRICH->GetPid(k));
      	  event.richbeta.push_back(hitRICH->GetBeta(k));
	}
      }

      //PID1
      const s_HodoRHitContainer &contPID1 =hit->GetsPID1RHC();
      int nhPID1 = contPID1.size();
      event.pid1nhits = nhPID1;
      for(int j=0; j<nhPID1; j++){
       	s_HodoRawHit *hitPID1=contPID1[j];
	event.pid1layer.push_back(hitPID1->LayerId());
      	event.pid1seg.push_back(hitPID1->SegmentId());

      	int nt = hitPID1->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.pid1time.push_back(hitPID1->GetTime(k));
      	  event.pid1edep.push_back(hitPID1->GetEdep(k));
      	  event.pid1path.push_back(hitPID1->GetPath(k));
      	  event.pid1p.push_back(hitPID1->GetMom(k));
      	  event.pid1posx.push_back(hitPID1->GetPosX(k));
      	  event.pid1posy.push_back(hitPID1->GetPosY(k));
      	  event.pid1pid.push_back(hitPID1->GetPid(k));
      	  event.pid1beta.push_back(hitPID1->GetBeta(k));
	}
      }

      //PID2
      const s_HodoRHitContainer &contPID2 =hit->GetsPID2RHC();
      int nhPID2 = contPID2.size();
      event.pid2nhits = nhPID2;
      for(int j=0; j<nhPID2; j++){
       	s_HodoRawHit *hitPID2=contPID2[j];
	event.pid2layer.push_back(hitPID2->LayerId());
      	event.pid2seg.push_back(hitPID2->SegmentId());

      	int nt = hitPID2->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.pid2time.push_back(hitPID2->GetTime(k));
      	  event.pid2edep.push_back(hitPID2->GetEdep(k));
      	  event.pid2path.push_back(hitPID2->GetPath(k));
      	  event.pid2p.push_back(hitPID2->GetMom(k));
      	  event.pid2posx.push_back(hitPID2->GetPosX(k));
      	  event.pid2posy.push_back(hitPID2->GetPosY(k));
      	  event.pid2pid.push_back(hitPID2->GetPid(k));
      	  event.pid2beta.push_back(hitPID2->GetBeta(k));
	}
      }

      //MF
      const s_HodoRHitContainer &contMF =hit->GetsMFRHC();
      int nhMF = contMF.size();
      event.mfnhits = nhMF;
      for(int j=0; j<nhMF; j++){
       	s_HodoRawHit *hitMF=contMF[j];
	event.mflayer.push_back(hitMF->LayerId());
      	event.mfseg.push_back(hitMF->SegmentId());

      	int nt = hitMF->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.mftime.push_back(hitMF->GetTime(k));
      	  event.mfedep.push_back(hitMF->GetEdep(k));
      	  event.mfpath.push_back(hitMF->GetPath(k));
      	  event.mfp.push_back(hitMF->GetMom(k));
      	  event.mfposx.push_back(hitMF->GetPosX(k));
      	  event.mfposy.push_back(hitMF->GetPosY(k));
      	  event.mfpid.push_back(hitMF->GetPid(k));
      	  event.mfbeta.push_back(hitMF->GetBeta(k));
	}
      }

      //VD
      const s_HodoRHitContainer &contVD =hit->GetsVDRHC();
      int nhVD = contVD.size();
      event.vdnhits = nhVD;
      for(int j=0; j<nhVD; j++){
       	s_HodoRawHit *hitVD=contVD[j];
	event.vdlayer.push_back(hitVD->LayerId());
      	event.vdseg.push_back(hitVD->SegmentId());

      	int nt = hitVD->GetSize();
	for( int k=0; k<nt; k++ ){
      	  event.vdtime.push_back(hitVD->GetTime(k));
      	  event.vdedep.push_back(hitVD->GetEdep(k));
      	  event.vdpath.push_back(hitVD->GetPath(k));
      	  event.vdp.push_back(hitVD->GetMom(k));
      	  event.vdposx.push_back(hitVD->GetPosX(k));
      	  event.vdposy.push_back(hitVD->GetPosY(k));
      	  event.vdpid.push_back(hitVD->GetPid(k));
      	  event.vdbeta.push_back(hitVD->GetBeta(k));
	}
      }

      //SFT
      for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      	const s_TrRHitContainer &contSFT = hit->GetsSFTRHC(layer);
      	int nhSFT=contSFT.size();
      	event.sftnhits = nhSFT;
      	for( int j=0; j<nhSFT; ++j ){
      	  s_TrRawHit *hitSFT=contSFT[j];
      	  event.sftlayer.push_back(hitSFT->LayerId());

	  //std::cout<< "sftlayer= " << hitSFT->LayerId() << std::endl;

      	  int nt = hitSFT->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.sftposx.push_back(hitSFT->GetPosX(k));
      	    event.sftposy.push_back(hitSFT->GetPosY(k));
      	    event.sftdl.push_back(hitSFT->GetDL(k));
      	  }
      	}
      }

      //AFT
      for( int layer=1; layer<=NumOfLayersAFT; ++layer ){
      	const s_TrRHitContainer &contAFT = hit->GetsAFTRHC(layer);
      	int nhAFT=contAFT.size();
      	event.aftnhits = nhAFT;
      	for( int j=0; j<nhAFT; ++j ){
      	  s_TrRawHit *hitAFT=contAFT[j];
      	  event.aftlayer.push_back(hitAFT->LayerId());

      	  int nt = hitAFT->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.aftposx.push_back(hitAFT->GetPosX(k));
      	    event.aftposy.push_back(hitAFT->GetPosY(k));
      	    event.aftdl.push_back(hitAFT->GetDL(k));
      	  }
      	}
      }

      //IT1
      for( int layer=1; layer<=NumOfLayersIT1; ++layer ){
      	const s_TrRHitContainer &contIT1 = hit->GetsIT1RHC(layer);
      	int nhIT1=contIT1.size();
      	event.it1nhits = nhIT1;
      	for( int j=0; j<nhIT1; ++j ){
      	  s_TrRawHit *hitIT1=contIT1[j];
      	  event.it1layer.push_back(hitIT1->LayerId());

	  //std::cout<< "it1layer= " << hitIT1->LayerId() << std::endl;

      	  int nt = hitIT1->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.it1posx.push_back(hitIT1->GetPosX(k));
      	    event.it1posy.push_back(hitIT1->GetPosY(k));
      	    event.it1dl.push_back(hitIT1->GetDL(k));
      	  }
      	}
      }

      //IT2
      for( int layer=1; layer<=NumOfLayersIT2; ++layer ){
      	const s_TrRHitContainer &contIT2 = hit->GetsIT2RHC(layer);
      	int nhIT2=contIT2.size();
      	event.it2nhits = nhIT2;
      	for( int j=0; j<nhIT2; ++j ){
      	  s_TrRawHit *hitIT2=contIT2[j];
      	  event.it2layer.push_back(hitIT2->LayerId());

      	  int nt = hitIT2->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.it2posx.push_back(hitIT2->GetPosX(k));
      	    event.it2posy.push_back(hitIT2->GetPosY(k));
      	    event.it2dl.push_back(hitIT2->GetDL(k));
      	  }
      	}
      }

      //ST1
      for( int layer=1; layer<=NumOfLayersST1; ++layer ){
      	const s_TrRHitContainer &contST1 = hit->GetsST1RHC(layer);
      	int nhST1=contST1.size();
      	event.st1nhits = nhST1;
      	for( int j=0; j<nhST1; ++j ){
      	  s_TrRawHit *hitST1=contST1[j];
      	  event.st1layer.push_back(hitST1->LayerId());

	  //std::cout<< "st1layer= " << hitST1->LayerId() << std::endl;

      	  int nt = hitST1->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.st1posx.push_back(hitST1->GetPosX(k));
      	    event.st1posy.push_back(hitST1->GetPosY(k));
      	    event.st1dl.push_back(hitST1->GetDL(k));
      	  }
      	}
      }

      //ST2
      for( int layer=1; layer<=NumOfLayersST2; ++layer ){
      	const s_TrRHitContainer &contST2 = hit->GetsST2RHC(layer);
      	int nhST2=contST2.size();
      	event.st2nhits = nhST2;
      	for( int j=0; j<nhST2; ++j ){
      	  s_TrRawHit *hitST2=contST2[j];
      	  event.st2layer.push_back(hitST2->LayerId());

	  //std::cout<< "st2layer= " << hitST2->LayerId() << std::endl;

      	  int nt = hitST2->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.st2posx.push_back(hitST2->GetPosX(k));
      	    event.st2posy.push_back(hitST2->GetPosY(k));
      	    event.st2dl.push_back(hitST2->GetDL(k));
      	  }
      	}
      }
    }
  }


  tree->Fill();

  return true;
}

void EventSMonitor::InitializeEvent( void )
{
  //PriInfo
  event.priposx = -9999.0;
  event.priposy = -9999.0;
  event.priposz = -9999.0;
  event.pbeam = -9999.0;
  event.ubeam = -9999.0;
  event.vbeam = -9999.0;
  event.abmom = -9999.0;
  event.prim1 = -9999.0;
  event.prip1 = -9999.0;
  event.pritheta1 = -9999.0;
  event.priphi1 = -9999.0;
  event.prithetacm1 = -9999.0;
  event.priphicm1 = -9999.0;
  event.prim2 = -9999.0;
  event.prip2 = -9999.0;
  event.pritheta2 = -9999.0;
  event.priphi2 = -9999.0;
  event.prithetacm2 = -9999.0;
  event.priphicm2 = -9999.0;

  //Generated Beam
  event.gbnhits = -1;
  event.gbid.clear();
  event.gbtype.clear();
  event.gbp.clear();
  event.gbpx.clear();
  event.gbpy.clear();
  event.gbpz.clear();
  event.gbpid.clear();
  event.gbvx.clear();
  event.gbvy.clear();

  //Generated Scattered
  event.gsnhits = -1;
  event.gsid.clear();
  event.gstype.clear();
  event.gsp.clear();
  event.gspx.clear();
  event.gspy.clear();
  event.gspz.clear();
  event.gspid.clear();
  event.gsvx.clear();
  event.gsvy.clear();

  //BFT
  event.bftnhits = -1;
  event.bftlayer.clear();
  event.bftposx.clear();
  event.bftposy.clear();
  event.bftdl.clear();

  //SFT
  event.sftnhits = -1;
  event.sftlayer.clear();
  event.sftposx.clear();
  event.sftposy.clear();
  event.sftdl.clear();

  //AFT
  event.aftnhits = -1;
  event.aftlayer.clear();
  event.aftposx.clear();
  event.aftposy.clear();
  event.aftdl.clear();

  //IT1
  event.it1nhits = -1;
  event.it1layer.clear();
  event.it1posx.clear();
  event.it1posy.clear();
  event.it1dl.clear();

  //IT2
  event.it2nhits = -1;
  event.it2layer.clear();
  event.it2posx.clear();
  event.it2posy.clear();
  event.it2dl.clear();

  //ST1
  event.st1nhits = -1;
  event.st1layer.clear();
  event.st1posx.clear();
  event.st1posy.clear();
  event.st1dl.clear();

  //ST2
  event.st2nhits = -1;
  event.st2layer.clear();
  event.st2posx.clear();
  event.st2posy.clear();
  event.st2dl.clear();

  //T0
  event.t0nhits = -1;
  event.t0layer.clear();
  event.t0seg.clear();
  event.t0time.clear();
  event.t0edep.clear();
  event.t0path.clear();
  event.t0p.clear();
  event.t0posx.clear();
  event.t0posy.clear();
  event.t0pid.clear();
  event.t0beta.clear();

  //TOF
  event.tofnhits = -1;
  event.toflayer.clear();
  event.tofseg.clear();
  event.toftime.clear();
  event.tofedep.clear();
  event.tofpath.clear();
  event.tofp.clear();
  event.tofposx.clear();
  event.tofposy.clear();
  event.tofpid.clear();
  event.tofbeta.clear();

  //ITOF
  event.itofnhits = -1;
  event.itoflayer.clear();
  event.itofseg.clear();
  event.itoftime.clear();
  event.itofedep.clear();
  event.itofpath.clear();
  event.itofp.clear();
  event.itofposx.clear();
  event.itofposy.clear();
  event.itofpid.clear();
  event.itofbeta.clear();

  //PAD
  event.padnhits = -1;
  event.padlayer.clear();
  event.padseg.clear();
  event.padtime.clear();
  event.padedep.clear();
  event.padpath.clear();
  event.padp.clear();
  event.padposx.clear();
  event.padposy.clear();
  event.padpid.clear();
  event.padbeta.clear();

  //RICH
  event.richnhits = -1;
  event.richlayer.clear();
  event.richseg.clear();
  event.richtime.clear();
  event.richedep.clear();
  event.richpath.clear();
  event.richp.clear();
  event.richposx.clear();
  event.richposy.clear();
  event.richpid.clear();
  event.richbeta.clear();

  //PID1
  event.pid1nhits = -1;
  event.pid1layer.clear();
  event.pid1seg.clear();
  event.pid1time.clear();
  event.pid1edep.clear();
  event.pid1path.clear();
  event.pid1p.clear();
  event.pid1posx.clear();
  event.pid1posy.clear();
  event.pid1pid.clear();
  event.pid1beta.clear();

  //PID2
  event.pid2nhits = -1;
  event.pid2layer.clear();
  event.pid2seg.clear();
  event.pid2time.clear();
  event.pid2edep.clear();
  event.pid2path.clear();
  event.pid2p.clear();
  event.pid2posx.clear();
  event.pid2posy.clear();
  event.pid2pid.clear();
  event.pid2beta.clear();

  //MF
  event.mfnhits = -1;
  event.mflayer.clear();
  event.mfseg.clear();
  event.mftime.clear();
  event.mfedep.clear();
  event.mfpath.clear();
  event.mfp.clear();
  event.mfposx.clear();
  event.mfposy.clear();
  event.mfpid.clear();
  event.mfbeta.clear();

  //VD
  event.vdnhits = -1;
  event.vdlayer.clear();
  event.vdseg.clear();
  event.vdtime.clear();
  event.vdedep.clear();
  event.vdpath.clear();
  event.vdp.clear();
  event.vdposx.clear();
  event.vdposy.clear();
  event.vdpid.clear();
  event.vdbeta.clear();
}


bool EventSMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSMonitor;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //PriInfo  
  tree->Branch("priposx", &event.priposx);
  tree->Branch("priposy", &event.priposy);
  tree->Branch("priposz", &event.priposz);
  tree->Branch("pbeam", &event.pbeam);
  tree->Branch("ubeam", &event.ubeam);
  tree->Branch("vbeam", &event.vbeam);
  tree->Branch("abmom", &event.abmom);
  tree->Branch("prim1", &event.prim1);
  tree->Branch("prip1", &event.prip1);
  tree->Branch("pritheta1", &event.pritheta1);
  tree->Branch("priphi1", &event.priphi1);
  tree->Branch("prithetacm1", &event.prithetacm1);
  tree->Branch("priphicm1", &event.priphicm1);
  tree->Branch("prim2", &event.prim2);
  tree->Branch("prip2", &event.prip2);
  tree->Branch("pritheta2", &event.pritheta2);
  tree->Branch("priphi2", &event.priphi2);
  tree->Branch("prithetacm2", &event.prithetacm2);
  tree->Branch("priphicm2", &event.priphicm2);

  //Generated Beam
  tree->Branch("gbnhits", &event.gbnhits);
  tree->Branch("gbid", &event.gbid);
  tree->Branch("gbtype", &event.gbtype);
  tree->Branch("gbp", &event.gbp);
  tree->Branch("gbpx", &event.gbpx);
  tree->Branch("gbpy", &event.gbpy);
  tree->Branch("gbpz", &event.gbpz);
  tree->Branch("gbpid", &event.gbpid);
  tree->Branch("gbvx", &event.gbvx);
  tree->Branch("gbvy", &event.gbvy);

  //Generated Scattered
  tree->Branch("gsnhits", &event.gsnhits);
  tree->Branch("gsid", &event.gsid);
  tree->Branch("gstype", &event.gstype);
  tree->Branch("gsp", &event.gsp);
  tree->Branch("gspx", &event.gspx);
  tree->Branch("gspy", &event.gspy);
  tree->Branch("gspz", &event.gspz);
  tree->Branch("gspid", &event.gspid);
  tree->Branch("gsvx", &event.gsvx);
  tree->Branch("gsvy", &event.gsvy);

  //BFT
  tree->Branch("bftnhits", &event.bftnhits);
  tree->Branch("bftlayer", &event.bftlayer);
  tree->Branch("bftposx", &event.bftposx);
  tree->Branch("bftposy", &event.bftposy);
  tree->Branch("bftdl", &event.bftdl);

  //SFT
  tree->Branch("sftnhits", &event.sftnhits);
  tree->Branch("sftlayer", &event.sftlayer);
  tree->Branch("sftposx", &event.sftposx);
  tree->Branch("sftposy", &event.sftposy);
  tree->Branch("sftdl", &event.sftdl);

  //AFT
  tree->Branch("aftnhits", &event.aftnhits);
  tree->Branch("aftlayer", &event.aftlayer);
  tree->Branch("aftposx", &event.aftposx);
  tree->Branch("aftposy", &event.aftposy);
  tree->Branch("aftdl", &event.aftdl);

  //IT1
  tree->Branch("it1nhits", &event.it1nhits);
  tree->Branch("it1layer", &event.it1layer);
  tree->Branch("it1posx", &event.it1posx);
  tree->Branch("it1posy", &event.it1posy);
  tree->Branch("it1dl", &event.it1dl);

  //IT2
  tree->Branch("it2nhits", &event.it2nhits);
  tree->Branch("it2layer", &event.it2layer);
  tree->Branch("it2posx", &event.it2posx);
  tree->Branch("it2posy", &event.it2posy);
  tree->Branch("it2dl", &event.it2dl);

  //ST1
  tree->Branch("st1nhits", &event.st1nhits);
  tree->Branch("st1layer", &event.st1layer);
  tree->Branch("st1posx", &event.st1posx);
  tree->Branch("st1posy", &event.st1posy);
  tree->Branch("st1dl", &event.st1dl);

  //ST2
  tree->Branch("st2nhits", &event.st2nhits);
  tree->Branch("st2layer", &event.st2layer);
  tree->Branch("st2posx", &event.st2posx);
  tree->Branch("st2posy", &event.st2posy);
  tree->Branch("st2dl", &event.st2dl);

  //T0
  tree->Branch("t0nhits", &event.t0nhits);
  tree->Branch("t0layer", &event.t0layer);
  tree->Branch("t0seg", &event.t0seg);
  tree->Branch("t0time", &event.t0time);
  tree->Branch("t0edep", &event.t0edep);
  tree->Branch("t0path", &event.t0path);
  tree->Branch("t0p", &event.t0p);
  tree->Branch("t0posx", &event.t0posx);
  tree->Branch("t0posy", &event.t0posy);
  tree->Branch("t0pid", &event.t0pid);
  tree->Branch("t0beta", &event.t0beta);
  
  //TOF
  tree->Branch("tofnhits", &event.tofnhits);
  tree->Branch("toflayer", &event.toflayer);
  tree->Branch("tofseg", &event.tofseg);
  tree->Branch("toftime", &event.toftime);
  tree->Branch("tofedep", &event.tofedep);
  tree->Branch("tofpath", &event.tofpath);
  tree->Branch("tofp", &event.tofp);
  tree->Branch("tofposx", &event.tofposx);
  tree->Branch("tofposy", &event.tofposy);
  tree->Branch("tofpid", &event.tofpid);
  tree->Branch("tofbeta", &event.tofbeta);
  
  //ITOF
  tree->Branch("itofnhits", &event.itofnhits);
  tree->Branch("itoflayer", &event.itoflayer);
  tree->Branch("itofseg", &event.itofseg);
  tree->Branch("itoftime", &event.itoftime);
  tree->Branch("itofedep", &event.itofedep);
  tree->Branch("itofpath", &event.itofpath);
  tree->Branch("itofp", &event.itofp);
  tree->Branch("itofposx", &event.itofposx);
  tree->Branch("itofposy", &event.itofposy);
  tree->Branch("itofpid", &event.itofpid);
  tree->Branch("itofbeta", &event.itofbeta);
  
  //PAD
  tree->Branch("padnhits", &event.padnhits);
  tree->Branch("padlayer", &event.padlayer);
  tree->Branch("padseg", &event.padseg);
  tree->Branch("padtime", &event.padtime);
  tree->Branch("padedep", &event.padedep);
  tree->Branch("padpath", &event.padpath);
  tree->Branch("padp", &event.padp);
  tree->Branch("padposx", &event.padposx);
  tree->Branch("padposy", &event.padposy);
  tree->Branch("padpid", &event.padpid);
  tree->Branch("padbeta", &event.padbeta);
  
  //RICH
  tree->Branch("richnhits", &event.richnhits);
  tree->Branch("richlayer", &event.richlayer);
  tree->Branch("richseg", &event.richseg);
  tree->Branch("richtime", &event.richtime);
  tree->Branch("richedep", &event.richedep);
  tree->Branch("richpath", &event.richpath);
  tree->Branch("richp", &event.richp);
  tree->Branch("richposx", &event.richposx);
  tree->Branch("richposy", &event.richposy);
  tree->Branch("richpid", &event.richpid);
  tree->Branch("richbeta", &event.richbeta);
  
  //PID1
  tree->Branch("pid1nhits", &event.pid1nhits);
  tree->Branch("pid1layer", &event.pid1layer);
  tree->Branch("pid1seg", &event.pid1seg);
  tree->Branch("pid1time", &event.pid1time);
  tree->Branch("pid1edep", &event.pid1edep);
  tree->Branch("pid1path", &event.pid1path);
  tree->Branch("pid1p", &event.pid1p);
  tree->Branch("pid1posx", &event.pid1posx);
  tree->Branch("pid1posy", &event.pid1posy);
  tree->Branch("pid1pid", &event.pid1pid);
  tree->Branch("pid1beta", &event.pid1beta);
  
  //PID2
  tree->Branch("pid2nhits", &event.pid2nhits);
  tree->Branch("pid2layer", &event.pid2layer);
  tree->Branch("pid2seg", &event.pid2seg);
  tree->Branch("pid2time", &event.pid2time);
  tree->Branch("pid2edep", &event.pid2edep);
  tree->Branch("pid2path", &event.pid2path);
  tree->Branch("pid2p", &event.pid2p);
  tree->Branch("pid2posx", &event.pid2posx);
  tree->Branch("pid2posy", &event.pid2posy);
  tree->Branch("pid2pid", &event.pid2pid);
  tree->Branch("pid2beta", &event.pid2beta);
  
  //MF
  tree->Branch("mfnhits", &event.mfnhits);
  tree->Branch("mflayer", &event.mflayer);
  tree->Branch("mfseg", &event.mfseg);
  tree->Branch("mftime", &event.mftime);
  tree->Branch("mfedep", &event.mfedep);
  tree->Branch("mfpath", &event.mfpath);
  tree->Branch("mfp", &event.mfp);
  tree->Branch("mfposx", &event.mfposx);
  tree->Branch("mfposy", &event.mfposy);
  tree->Branch("mfpid", &event.mfpid);
  tree->Branch("mfbeta", &event.mfbeta);
  
  //VD
  tree->Branch("vdnhits", &event.vdnhits);
  tree->Branch("vdlayer", &event.vdlayer);
  tree->Branch("vdseg", &event.vdseg);
  tree->Branch("vdtime", &event.vdtime);
  tree->Branch("vdedep", &event.vdedep);
  tree->Branch("vdpath", &event.vdpath);
  tree->Branch("vdp", &event.vdp);
  tree->Branch("vdposx", &event.vdposx);
  tree->Branch("vdposy", &event.vdposy);
  tree->Branch("vdpid", &event.vdpid);
  tree->Branch("vdbeta", &event.vdbeta);

  return true;
}
