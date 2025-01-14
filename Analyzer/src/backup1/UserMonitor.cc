/*
  UserMonitor.cc

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

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventMonitor();
  ~EventMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventMonitor::EventMonitor()
  : VEvent(),
    rawData(0)
{
}

EventMonitor::~EventMonitor()
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

bool EventMonitor::ProcessingBegin()
{
 return true;
}

bool EventMonitor::ProcessingNormal( std::ifstream &In )
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

  //BFT
  {
    for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetBFTRHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.bftnhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.bftlayer.push_back(hit->LayerId());
  	  event.bftposx.push_back(hit->GetPosX(j));
  	  event.bftposy.push_back(hit->GetPosY(j));
  	  event.bftdl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //SFT
  {
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetSFTRHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.sftnhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.sftlayer.push_back(hit->LayerId());
  	  event.sftposx.push_back(hit->GetPosX(j));
  	  event.sftposy.push_back(hit->GetPosY(j));
  	  event.sftdl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //AFT
  {
    for( int layer=1; layer<=NumOfLayersAFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetAFTRHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.aftnhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.aftlayer.push_back(hit->LayerId());
  	  event.aftposx.push_back(hit->GetPosX(j));
  	  event.aftposy.push_back(hit->GetPosY(j));
  	  event.aftdl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //IT1
  {
    for( int layer=1; layer<=NumOfLayersIT1; ++layer ){
      const TrRHitContainer &cont =rawData->GetIT1RHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.it1nhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.it1layer.push_back(hit->LayerId());
  	  event.it1posx.push_back(hit->GetPosX(j));
  	  event.it1posy.push_back(hit->GetPosY(j));
  	  event.it1dl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //IT2
  {
    for( int layer=1; layer<=NumOfLayersIT2; ++layer ){
      const TrRHitContainer &cont =rawData->GetIT2RHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.it2nhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.it2layer.push_back(hit->LayerId());
  	  event.it2posx.push_back(hit->GetPosX(j));
  	  event.it2posy.push_back(hit->GetPosY(j));
  	  event.it2dl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //ST1
  {
    for( int layer=1; layer<=NumOfLayersST1; ++layer ){
      const TrRHitContainer &cont =rawData->GetST1RHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.st1nhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.st1layer.push_back(hit->LayerId());
  	  event.st1posx.push_back(hit->GetPosX(j));
  	  event.st1posy.push_back(hit->GetPosY(j));
  	  event.st1dl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //ST2
  {
    for( int layer=1; layer<=NumOfLayersST2; ++layer ){
      const TrRHitContainer &cont =rawData->GetST2RHC(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.st2nhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.st2layer.push_back(hit->LayerId());
  	  event.st2posx.push_back(hit->GetPosX(j));
  	  event.st2posy.push_back(hit->GetPosY(j));
  	  event.st2dl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //T0
  {
    const HodoRHitContainer &cont=rawData->GetT0RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.t0nhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.t0layer.push_back(layer);
  	event.t0seg.push_back(seg);
  	event.t0time.push_back(hit->GetTime(j));
  	event.t0edep.push_back(hit->GetEdep(j));
  	event.t0path.push_back(hit->GetPath(j));
  	event.t0p.push_back(hit->GetMom(j));
  	event.t0posx.push_back(hit->GetPosX(j));
  	event.t0posy.push_back(hit->GetPosY(j));
  	event.t0pid.push_back(hit->GetPid(j));
  	event.t0beta.push_back(hit->GetBeta(j));
      }
    }
  }

  //TOF
  {
    const HodoRHitContainer &cont=rawData->GetTOFRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.tofnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.toflayer.push_back(layer);
  	event.tofseg.push_back(seg);
  	event.toftime.push_back(hit->GetTime(j));
  	event.tofedep.push_back(hit->GetEdep(j));
  	event.tofpath.push_back(hit->GetPath(j));
  	event.tofp.push_back(hit->GetMom(j));
  	event.tofposx.push_back(hit->GetPosX(j));
  	event.tofposy.push_back(hit->GetPosY(j));
  	event.tofpid.push_back(hit->GetPid(j));
  	event.tofbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  //ITOF
  {
    const HodoRHitContainer &cont=rawData->GetITOFRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.itofnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.itoflayer.push_back(layer);
  	event.itofseg.push_back(seg);
  	event.itoftime.push_back(hit->GetTime(j));
  	event.itofedep.push_back(hit->GetEdep(j));
  	event.itofpath.push_back(hit->GetPath(j));
  	event.itofp.push_back(hit->GetMom(j));
  	event.itofposx.push_back(hit->GetPosX(j));
  	event.itofposy.push_back(hit->GetPosY(j));
  	event.itofpid.push_back(hit->GetPid(j));
  	event.itofbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  //PAD
  {
    const HodoRHitContainer &cont=rawData->GetPADRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.padnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.padlayer.push_back(layer);
  	event.padseg.push_back(seg);
  	event.padtime.push_back(hit->GetTime(j));
  	event.padedep.push_back(hit->GetEdep(j));
  	event.padpath.push_back(hit->GetPath(j));
  	event.padp.push_back(hit->GetMom(j));
  	event.padposx.push_back(hit->GetPosX(j));
  	event.padposy.push_back(hit->GetPosY(j));
  	event.padpid.push_back(hit->GetPid(j));
  	event.padbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  //RICH
  {
    const HodoRHitContainer &cont=rawData->GetRICHRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.richnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.richlayer.push_back(layer);
  	event.richseg.push_back(seg);
  	event.richtime.push_back(hit->GetTime(j));
  	event.richedep.push_back(hit->GetEdep(j));
  	event.richpath.push_back(hit->GetPath(j));
  	event.richp.push_back(hit->GetMom(j));
  	event.richposx.push_back(hit->GetPosX(j));
  	event.richposy.push_back(hit->GetPosY(j));
  	event.richpid.push_back(hit->GetPid(j));
  	event.richbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  //PID1
  {
    const HodoRHitContainer &cont=rawData->GetPID1RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.pid1nhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.pid1layer.push_back(layer);
  	event.pid1seg.push_back(seg);
  	event.pid1time.push_back(hit->GetTime(j));
  	event.pid1edep.push_back(hit->GetEdep(j));
  	event.pid1path.push_back(hit->GetPath(j));
  	event.pid1p.push_back(hit->GetMom(j));
  	event.pid1posx.push_back(hit->GetPosX(j));
  	event.pid1posy.push_back(hit->GetPosY(j));
  	event.pid1pid.push_back(hit->GetPid(j));
  	event.pid1beta.push_back(hit->GetBeta(j));
      }
    }
  }

  //PID2
  {
    const HodoRHitContainer &cont=rawData->GetPID2RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.pid2nhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.pid2layer.push_back(layer);
  	event.pid2seg.push_back(seg);
  	event.pid2time.push_back(hit->GetTime(j));
  	event.pid2edep.push_back(hit->GetEdep(j));
  	event.pid2path.push_back(hit->GetPath(j));
  	event.pid2p.push_back(hit->GetMom(j));
  	event.pid2posx.push_back(hit->GetPosX(j));
  	event.pid2posy.push_back(hit->GetPosY(j));
  	event.pid2pid.push_back(hit->GetPid(j));
  	event.pid2beta.push_back(hit->GetBeta(j));
      }
    }
  }

  //MF
  {
    const HodoRHitContainer &cont=rawData->GetMFRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.mfnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.mflayer.push_back(layer);
  	event.mfseg.push_back(seg);
  	event.mftime.push_back(hit->GetTime(j));
  	event.mfedep.push_back(hit->GetEdep(j));
  	event.mfpath.push_back(hit->GetPath(j));
  	event.mfp.push_back(hit->GetMom(j));
  	event.mfposx.push_back(hit->GetPosX(j));
  	event.mfposy.push_back(hit->GetPosY(j));
  	event.mfpid.push_back(hit->GetPid(j));
  	event.mfbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  //VD
  {
    const HodoRHitContainer &cont=rawData->GetVDRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.vdnhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.vdlayer.push_back(layer);
  	event.vdseg.push_back(seg);
  	event.vdtime.push_back(hit->GetTime(j));
  	event.vdedep.push_back(hit->GetEdep(j));
  	event.vdpath.push_back(hit->GetPath(j));
  	event.vdp.push_back(hit->GetMom(j));
  	event.vdposx.push_back(hit->GetPosX(j));
  	event.vdposy.push_back(hit->GetPosY(j));
  	event.vdpid.push_back(hit->GetPid(j));
  	event.vdbeta.push_back(hit->GetBeta(j));
      }
    }
  }

  tree->Fill();

  return true;
}

void EventMonitor::InitializeEvent( void )
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


bool EventMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor;
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
