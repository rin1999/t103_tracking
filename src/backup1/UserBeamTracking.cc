/*
  UserBeamTracking.cc

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
#include "SpecLib.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBeamTracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  //  HodoAnalyzer *HodoAna;

public:
  EventBeamTracking();
  ~EventBeamTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventBeamTracking::EventBeamTracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer())
    //HodoAna(new HodoAnalyzer())
{
}

EventBeamTracking::~EventBeamTracking()
{
  //  if (HodoAna)   delete HodoAna;
  if (TrAna)   delete TrAna;
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

  //BFT
  int    bftnhits;
  std::vector<int> bftlayer;
  std::vector<double> bftposx, bftposy;
  std::vector<double> bftdl;

  //T0
  int    t0nhits;
  std::vector<int>    t0layer, t0seg;
  std::vector<double> t0time, t0edep;
  std::vector<double> t0path, t0p;
  std::vector<double> t0posx, t0posy;
  std::vector<int>    t0pid;
  std::vector<double> t0beta;

  //Local tracking
  int    ntbl;
  std::vector<int>    idbl, typebl;
  std::vector<int>    layerbl;
  std::vector<double> chisqrbl;
  std::vector<double> x0bl, y0bl;
  std::vector<double> u0bl, v0bl;
  std::vector<double> posbl, resbl;

  //Beam Tracking
  int    ntb;
  std::vector<int>    idb, typeb;
  std::vector<int>    layerb;
  std::vector<double> chisqrb;
  std::vector<double> x0b, y0b;
  std::vector<double> u0b, v0b;
  std::vector<double> posb, resb;

};
static Event event;

bool EventBeamTracking::ProcessingBegin()
{
 return true;
}

bool EventBeamTracking::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits(In) ) return false;
  //std::cout << "***" << std::endl;

  //**************************************************************************
  //******************RawData

  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  double X0, Y0, Z0;
  //PrimaryInfo
  {
    const PrimInfoContainer &cont=rawData->GetPrimHC();  
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      PrimInfo *hit=cont[i];
      event.priposx = hit->GetVertX();
      event.priposy = hit->GetVertY();
      event.priposz = hit->GetVertZ();
      X0=hit->GetVertX();
      Y0=hit->GetVertY();
      Z0=hit->GetVertZ();
      
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
  ThreeVector IniVert(1400.0-Z0, 0.0, 0.0);

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
      	s_HodoRawHit *hitT0=contT0[j];
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
      	int nhBFT=contBFT.size();
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

      //////////////////////////Tracking
      event.idbl.push_back(hit->TrackId());
      event.typebl.push_back(hit->TrackType());
      
      TrAna->DecodesBRawHits( hit, hit->TrackType() );

      //BFT
      TrAna->TrackSearchBFTT();
      int nt=TrAna->GetNtracksBFTT();
      event.ntbl=nt;
      for( int it=0; it<nt; ++it ){
	TrLocalTrack *tp=TrAna->GetTrackBFTT(it);
	int nh=tp->GetNHit();
	double chisqr=tp->GetChiSquare();
	double x0=tp->GetX0(), y0=tp->GetY0();
	double u0=tp->GetU0(), v0=tp->GetV0();
	double xtgt=tp->GetX( Z0 ), ytgt=tp->GetY( Z0 );
	double utgt=u0, vtgt=v0;

	event.chisqrbl.push_back(chisqr);
	event.x0bl.push_back(xtgt);
	event.y0bl.push_back(ytgt);
	event.u0bl.push_back(utgt);
	event.v0bl.push_back(vtgt); 
	
	for( int ih=0; ih<nh; ++ih ){
	  TrLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer(); 
	  event.layerbl.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.posbl.push_back(pos);
	  event.resbl.push_back(res);
	}
      }

      TrAna->TrackSearchBeam( IniVert );
      int ntBeam=TrAna->GetNTracksBeam();
      event.ntb = ntBeam;
      for( int it=0; it<ntBeam; ++it ){
	BeamTrack *tp=TrAna->GetBeamTrack(it);
	if(!tp) continue;
	int nh=tp->GetNHits();
	double chisqr=tp->chisqr();  
	ThreeVector Ppos=tp->PrimaryPosition();
	ThreeVector Pmom=tp->PrimaryMomentum();
	double x=Ppos.x(), y=Ppos.y();
	double p=Pmom.mag();
	double u=Pmom.x()/p, v=Pmom.y()/p;
	
	// std::cout<< "V1= " << ThreeVector(X0, Y0, 0.0) << std::endl;
	// std::cout<< "V2= " << Ppos << std::endl;
	//std::cout<< Pmom << std::endl;

	event.chisqrb.push_back(chisqr);
	event.x0b.push_back(x);
	event.y0b.push_back(y);
	event.u0b.push_back(u);
	event.v0b.push_back(v);   

	for( int j=0; j<nh; ++j ){
	  TrackHit *hit=tp->GetHit(j);
	  if(!hit) continue;
	  int layerId=hit->GetLayer();
	  event.layerb.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.posb.push_back(pos);
	  event.resb.push_back(res);
	}
      }
    }
  }


  tree->Fill();

  return true;
}

void EventBeamTracking::InitializeEvent( void )
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

  //BFT
  event.bftnhits = -1;
  event.bftlayer.clear();
  event.bftposx.clear();
  event.bftposy.clear();
  event.bftdl.clear();

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

  //Local Tracking
  event.idbl.clear();
  event.typebl.clear();
  event.ntbl = -1;
  event.layerbl.clear();
  event.chisqrbl.clear();
  event.x0bl.clear();
  event.y0bl.clear();
  event.u0bl.clear();
  event.v0bl.clear();
  event.posbl.clear();
  event.resbl.clear();

  //Beam Tracking
  event.idb.clear();
  event.typeb.clear();
  event.ntb = -1;
  event.layerb.clear();
  event.chisqrb.clear();
  event.x0b.clear();
  event.y0b.clear();
  event.u0b.clear();
  event.v0b.clear();
  event.posb.clear();
  event.resb.clear();
}

bool EventBeamTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBeamTracking;
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

  //BFT
  tree->Branch("bftnhits", &event.bftnhits);
  tree->Branch("bftlayer", &event.bftlayer);
  tree->Branch("bftposx", &event.bftposx);
  tree->Branch("bftposy", &event.bftposy);
  tree->Branch("bftdl", &event.bftdl);

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
  
  //Local Tracking
  tree->Branch("idbl", &event.idbl);
  tree->Branch("typebl", &event.typebl);
  tree->Branch("ntbl", &event.ntbl);
  tree->Branch("layerbl", &event.layerbl);
  tree->Branch("chisqrbl", &event.chisqrbl);
  tree->Branch("x0bl", &event.x0bl);
  tree->Branch("y0bl", &event.y0bl);
  tree->Branch("u0bl", &event.u0bl);
  tree->Branch("v0bl", &event.v0bl);
  tree->Branch("posbl", &event.posbl);
  tree->Branch("resbl", &event.resbl);

  //Beam Tracking
  tree->Branch("idb", &event.idb);
  tree->Branch("typeb", &event.typeb);
  tree->Branch("ntb", &event.ntb);
  tree->Branch("layerb", &event.layerb);
  tree->Branch("chisqrb", &event.chisqrb);
  tree->Branch("x0b", &event.x0b);
  tree->Branch("y0b", &event.y0b);
  tree->Branch("u0b", &event.u0b);
  tree->Branch("v0b", &event.v0b);
  tree->Branch("posb", &event.posb);
  tree->Branch("resb", &event.resb);

  return true;
}
