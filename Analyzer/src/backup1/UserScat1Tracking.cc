/*
  UserScat1Tracking.cc

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

#define check 1

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventScat1Tracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  HodoAnalyzer *HodoAna;

public:
  EventScat1Tracking();
  ~EventScat1Tracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventScat1Tracking::EventScat1Tracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer()),
    HodoAna(new HodoAnalyzer())
{
}

EventScat1Tracking::~EventScat1Tracking()
{
  if (HodoAna)   delete HodoAna;
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

  //Generated Scat
  int    gsnhits;
  std::vector<int>    gsid, gstype;
  std::vector<double> gsp, gspx, gspy, gspz;
  std::vector<int>    gspid;
  std::vector<double> gsvx, gsvy;

  //SFT
  std::vector<int> sftnhits;
  std::vector<int> sftlayer;
  std::vector<double> sftposx, sftposy;
  std::vector<double> sftdl;

  //IT1
  std::vector<int> it1nhits;
  std::vector<int> it1layer;
  std::vector<double> it1posx, it1posy;
  std::vector<double> it1dl;

  //ST1
  std::vector<int> st1nhits;
  std::vector<int> st1layer;
  std::vector<double> st1posx, st1posy;
  std::vector<double> st1dl;

  //ST2
  std::vector<int> st2nhits;
  std::vector<int> st2layer;
  std::vector<double> st2posx, st2posy;
  std::vector<double> st2dl;

  //TOF
  std::vector<int>    tofnhits;
  std::vector<int>    toflayer, tofseg;
  std::vector<double> toftime, tofedep;
  std::vector<double> tofpath, tofp;
  std::vector<double> tofposx, tofposy;
  std::vector<int>    tofpid;
  std::vector<double> tofbeta;

  //RICH
  std::vector<int>    richnhits;
  std::vector<int>    richlayer, richseg;
  std::vector<double> richtime, richedep;
  std::vector<double> richpath, richp;
  std::vector<double> richposx, richposy;
  std::vector<int>    richpid;
  std::vector<double> richbeta;

  //PID1
  std::vector<int>    pid1nhits;
  std::vector<int>    pid1layer, pid1seg;
  std::vector<double> pid1time, pid1edep;
  std::vector<double> pid1path, pid1p;
  std::vector<double> pid1posx, pid1posy;
  std::vector<int>    pid1pid;
  std::vector<double> pid1beta;

  //Local tracking 
  std::vector<int>    ntslin1;
  std::vector<int>    layerslin1;
  std::vector<double> chisqrslin1;
  std::vector<double> x0slin1, y0slin1;
  std::vector<double> u0slin1, v0slin1;
  std::vector<double> posslin1, resslin1;

  std::vector<int>    ntslin2;
  std::vector<int>    layerslin2;
  std::vector<double> chisqrslin2;
  std::vector<double> x0slin2, y0slin2;
  std::vector<double> u0slin2, v0slin2;
  std::vector<double> posslin2, resslin2;

  std::vector<int>    ntslout1;
  std::vector<int>    layerslout1;
  std::vector<double> chisqrslout1;
  std::vector<double> x0slout1, y0slout1;
  std::vector<double> u0slout1, v0slout1;
  std::vector<double> posslout1, resslout1;

  std::vector<int>    ntslout2;
  std::vector<int>    layerslout2;
  std::vector<double> chisqrslout2;
  std::vector<double> x0slout2, y0slout2;
  std::vector<double> u0slout2, v0slout2;
  std::vector<double> posslout2, resslout2;

  //Scat Tracking
  std::vector<int>    ntsprein, ntspreout;
  std::vector<int>    nts;
  std::vector<int>    ids, types;
  std::vector<int>    layers;
  std::vector<double> chisqrs;
  std::vector<double> p, theta, phi;
  std::vector<double> dp;
  std::vector<double> path, m2;
  std::vector<double> x0s, y0s;
  std::vector<double> u0s, v0s;
  std::vector<double> poss, ress;
};
static Event event;

bool EventScat1Tracking::ProcessingBegin()
{
 return true;
}

bool EventScat1Tracking::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits(In) ) return false;

#if check
  std::cout << "*************************************************" << std::endl; 
#endif

 const TrGeomMan & geomMan=TrGeomMan::GetInstance();

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

  //Scat
  double IniP;
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

      //TOF
      const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
      int nhTOF = contTOF.size();
      event.tofnhits.push_back(nhTOF);
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

      //RICH
      const s_HodoRHitContainer &contRICH =hit->GetsRICHRHC();
      int nhRICH = contRICH.size();
      event.richnhits.push_back(nhRICH);
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
      event.pid1nhits.push_back(nhPID1);
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

      //SFT
      for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      	const s_TrRHitContainer &contSFT = hit->GetsSFTRHC(layer);
      	int nhSFT=contSFT.size();
      	event.sftnhits.push_back(nhSFT);
      	for( int j=0; j<nhSFT; ++j ){
      	  s_TrRawHit *hitSFT=contSFT[j];
      	  event.sftlayer.push_back(hitSFT->LayerId());

      	  int nt = hitSFT->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.sftposx.push_back(hitSFT->GetPosX(k));
      	    event.sftposy.push_back(hitSFT->GetPosY(k));
      	    event.sftdl.push_back(hitSFT->GetDL(k));
	
      	  }
      	}
      }

      //IT1
      for( int layer=1; layer<=NumOfLayersIT1; ++layer ){
      	const s_TrRHitContainer &contIT1 = hit->GetsIT1RHC(layer);
      	int nhIT1=contIT1.size();
      	event.it1nhits.push_back(nhIT1);
      	for( int j=0; j<nhIT1; ++j ){
      	  s_TrRawHit *hitIT1=contIT1[j];
      	  event.it1layer.push_back(hitIT1->LayerId());

      	  int nt = hitIT1->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.it1posx.push_back(hitIT1->GetPosX(k));
      	    event.it1posy.push_back(hitIT1->GetPosY(k));
      	    event.it1dl.push_back(hitIT1->GetDL(k));
      	  }
      	}
      }

      //ST1
      for( int layer=1; layer<=NumOfLayersST1; ++layer ){
      	const s_TrRHitContainer &contST1 = hit->GetsST1RHC(layer);
      	int nhST1=contST1.size();
      	event.st1nhits.push_back(nhST1);
      	for( int j=0; j<nhST1; ++j ){
      	  s_TrRawHit *hitST1=contST1[j];
      	  event.st1layer.push_back(hitST1->LayerId());

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
      	event.st2nhits.push_back(nhST2);
      	for( int j=0; j<nhST2; ++j ){
      	  s_TrRawHit *hitST2=contST2[j];
      	  event.st2layer.push_back(hitST2->LayerId());

      	  int nt = hitST2->GetSize();
      	  for( int k=0; k<nt; k++ ) {
      	    event.st2posx.push_back(hitST2->GetPosX(k));
      	    event.st2posy.push_back(hitST2->GetPosY(k));
      	    event.st2dl.push_back(hitST2->GetDL(k));
      	  }
      	}
      }

      //////////////////////////Tracking
      event.ids.push_back(hit->TrackId());
      event.types.push_back(hit->TrackType());
      IniP=hit->GetMom();

      if( hit->TrackType() == TrackTypeScat1 ){
	//std::cout<< "*****Before Tracking" << std::endl;
	TrAna->DecodesSRawHits( hit, hit->TrackType() );

	//SFT
	TrAna->TrackSearchSFTT();
	int ntin1=TrAna->GetNtracksSFTT();
	event.ntslin1.push_back(ntin1);
	for( int it=0; it<ntin1; ++it ){
	  TrLocalTrack *tp=TrAna->GetTrackSFTT(it);
	  int nh=tp->GetNHit();
	  double chisqr=tp->GetChiSquare();
	  double x0=tp->GetX0(), y0=tp->GetY0();
	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
	  int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	  double utgt=u0, vtgt=v0;
	  
	  event.chisqrslin1.push_back(chisqr);
	  event.x0slin1.push_back(xtgt);
	  event.y0slin1.push_back(ytgt);
	  event.u0slin1.push_back(utgt);
	  event.v0slin1.push_back(vtgt); 
	  
	  for( int ih=0; ih<nh; ++ih ){
	    TrLTrackHit *hit=tp->GetHit(ih);
	    int layerId=hit->GetLayer(); 
	    event.layerslin1.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.posslin1.push_back(pos);
	    event.resslin1.push_back(res);
	  }
	}
	
	//IT1
	TrAna->TrackSearchIT1T();
	int ntin2=TrAna->GetNtracksIT1T();
	event.ntslin2.push_back(ntin2);
	for( int it=0; it<ntin2; ++it ){
	  TrLocalTrack *tp=TrAna->GetTrackIT1T(it);
	  int nh=tp->GetNHit();
	  double chisqr=tp->GetChiSquare();
	  double x0=tp->GetX0(), y0=tp->GetY0();
	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
	  int Id = TrGeomMan::GetInstance().GetDetectorId("SpecVp2");
	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	  double utgt=u0, vtgt=v0;
	  
	  event.chisqrslin2.push_back(chisqr);
	  event.x0slin2.push_back(xtgt);
	  event.y0slin2.push_back(ytgt);
	  event.u0slin2.push_back(utgt);
	  event.v0slin2.push_back(vtgt); 
	  
	  for( int ih=0; ih<nh; ++ih ){
	    TrLTrackHit *hit=tp->GetHit(ih);
	    int layerId=hit->GetLayer(); 
	    event.layerslin2.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.posslin2.push_back(pos);
	    event.resslin2.push_back(res);
	  }
	}
	
	//ST1
	TrAna->TrackSearchST1T();
	int ntout1=TrAna->GetNtracksST1T();
	event.ntslout1.push_back(ntout1);
	for( int it=0; it<ntout1; ++it ){
	  TrLocalTrack *tp=TrAna->GetTrackST1T(it);
	  int nh=tp->GetNHit();
	  double chisqr=tp->GetChiSquare();
	  double x0=tp->GetX0(), y0=tp->GetY0();
	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
	  int Id = TrGeomMan::GetInstance().GetDetectorId("SpecVp3");
	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	  double utgt=u0, vtgt=v0;
	  
	  event.chisqrslout1.push_back(chisqr);
	  event.x0slout1.push_back(xtgt);
	  event.y0slout1.push_back(ytgt);
	  event.u0slout1.push_back(utgt);
	  event.v0slout1.push_back(vtgt); 
	  
	  for( int ih=0; ih<nh; ++ih ){
	    TrLTrackHit *hit=tp->GetHit(ih);
	    int layerId=hit->GetLayer(); 
	    event.layerslout1.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.posslout1.push_back(pos);
	    event.resslout1.push_back(res);
	  }
	}
	
	//ST2
	TrAna->TrackSearchST2T();
	int ntout2=TrAna->GetNtracksST2T();
	event.ntslout2.push_back(ntout2);
	for( int it=0; it<ntout2; ++it ){
	  TrLocalTrack *tp=TrAna->GetTrackST2T(it);
	  int nh=tp->GetNHit();
	  double chisqr=tp->GetChiSquare();
	  double x0=tp->GetX0(), y0=tp->GetY0();
	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
	  int Id = TrGeomMan::GetInstance().GetDetectorId("TOF");
	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	  double utgt=u0, vtgt=v0;
	  
	  event.chisqrslout2.push_back(chisqr);
	  event.x0slout2.push_back(xtgt);
	  event.y0slout2.push_back(ytgt);
	  event.u0slout2.push_back(utgt);
	  event.v0slout2.push_back(vtgt); 
	  
	  for( int ih=0; ih<nh; ++ih ){
	    TrLTrackHit *hit=tp->GetHit(ih);
	    int layerId=hit->GetLayer(); 
	    event.layerslout2.push_back(layerId);  
	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	    event.posslout2.push_back(pos);
	    event.resslout2.push_back(res);
	  }
	}
	
	//Scat1 tracking
      	TrAna->TrackSearchPreIn( IniP );
	TrAna->TrackSearchPreOut( IniP );

	int ntPreIn=TrAna->GetNTracksPreIn();
	int ntPreOut=TrAna->GetNTracksPreOut();
      	event.ntsprein.push_back(ntPreIn);
      	event.ntspreout.push_back(ntPreOut);

      	TrAna->TrackSearchScat1( IniP, IniVert );
      	int ntScat1=TrAna->GetNTracksScat1();
      	event.nts.push_back(ntScat1);
      	for( int it=0; it<ntScat1; ++it ){
      	  Scat1Track *tp=TrAna->GetScat1Track(it);
      	  if(!tp) continue;
      	  int nh=tp->GetNHits();
      	  double chisqr=tp->chisqr();  
      	  ThreeVector Ppos=tp->PrimaryPosition();
      	  ThreeVector Pmom=tp->PrimaryMomentum();
      	  double pathL=tp->PathLengthToTOF();
      	  double x=Ppos.x(), y=Ppos.y();
      	  double p=Pmom.mag();
      	  double pz=Pmom.z();
	  double u, v;
      	  double theta, phi;	  

	  if(pz>0){
	    u=Pmom.x()/pz, v=Pmom.y()/pz;
	    theta = Pmom.theta();	  
	    phi = Pmom.phi();	  
	  }
	  if(pz<0){
	    u=Pmom.x()/pz, v=Pmom.y()/pz;
	    theta = (-1.*Pmom).theta();	  
	    phi = (-1.*Pmom).phi();	  
	  }      

      	  // std::cout<< "V1= " << ThreeVector(X0, Y0, 0.0) << std::endl;
      	  // std::cout<< "V2= " << Ppos << std::endl;
	  
      	  event.chisqrs.push_back(chisqr);
      	  event.p.push_back(p);
	  event.theta.push_back(theta*Rad2Deg);
	  event.phi.push_back(phi*Rad2Deg);
      	  event.path.push_back(pathL);
      	  event.x0s.push_back(x);
      	  event.y0s.push_back(y);
      	  event.u0s.push_back(u);
      	  event.v0s.push_back(v);   

	  event.dp.push_back(p-IniP);

#if check
	  std::cout<< "***** p= " << p << ", dp= " << p-IniP <<std::endl;
#endif

	  //TOF
	  const s_HodoRHitContainer &contTOF =hit->GetsTOFRHC();
	  int nhTOF = contTOF.size();
	  for(int j=0; j<nhTOF; j++){
	    s_HodoRawHit *hitTOF=contTOF[j];
	    
	    int nt = hitTOF->GetSize();
	    for( int k=0; k<nt; k++ ){
	      double time = hitTOF->GetTime(k);
	      double path = hitTOF->GetPath(k);
	      double mom = hitTOF->GetMom(k);
	      
	      double C = 2.99792458E+8;
	      double T0, V0, B0, m2_0;
	      
	      V0 = (path/1000.)/(time*1.0E-9);
	      B0 = V0/C;
	      m2_0 = (mom/B0)*(mom/B0)*(1.-B0*B0);
	      
	      event.m2.push_back(m2_0); 
	    }
	  }
	  
	  for( int j=0; j<nh; ++j ){
	    TrackHit *hit=tp->GetHit(j);
	    if(!hit) continue;
	    int layerId=hit->GetLayer();
      	    event.layers.push_back(layerId);  
      	  }
	}
      }
    }
  }


  tree->Fill();

  return true;
}

void EventScat1Tracking::InitializeEvent( void )
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

  //Generated Scat
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

  //SFT
  event.sftnhits.clear();
  event.sftlayer.clear();
  event.sftposx.clear();
  event.sftposy.clear();
  event.sftdl.clear();

  //IT1
  event.it1nhits.clear();
  event.it1layer.clear();
  event.it1posx.clear();
  event.it1posy.clear();
  event.it1dl.clear();

  //ST1
  event.st1nhits.clear();
  event.st1layer.clear();
  event.st1posx.clear();
  event.st1posy.clear();
  event.st1dl.clear();

  //ST2
  event.st2nhits.clear();
  event.st2layer.clear();
  event.st2posx.clear();
  event.st2posy.clear();
  event.st2dl.clear();

  //TOF
  event.tofnhits.clear();
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

  //RICH
  event.richnhits.clear();
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
  event.pid1nhits.clear();
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

  //Local Tracking
  event.ntslin1.clear();
  event.layerslin1.clear();
  event.chisqrslin1.clear();
  event.x0slin1.clear();
  event.y0slin1.clear();
  event.u0slin1.clear();
  event.v0slin1.clear();
  event.posslin1.clear();
  event.resslin1.clear();

  event.ntslin2.clear();
  event.layerslin2.clear();
  event.chisqrslin2.clear();
  event.x0slin2.clear();
  event.y0slin2.clear();
  event.u0slin2.clear();
  event.v0slin2.clear();
  event.posslin2.clear();
  event.resslin2.clear();

  event.ntslout1.clear();
  event.layerslout1.clear();
  event.chisqrslout1.clear();
  event.x0slout1.clear();
  event.y0slout1.clear();
  event.u0slout1.clear();
  event.v0slout1.clear();
  event.posslout1.clear();
  event.resslout1.clear();

  event.ntslout2.clear();
  event.layerslout2.clear();
  event.chisqrslout2.clear();
  event.x0slout2.clear();
  event.y0slout2.clear();
  event.u0slout2.clear();
  event.v0slout2.clear();
  event.posslout2.clear();
  event.resslout2.clear();

  //Scat Tracking
  event.ids.clear();
  event.types.clear();
  event.ntsprein.clear();
  event.ntspreout.clear();
  event.nts.clear();
  event.layers.clear();
  event.chisqrs.clear();
  event.p.clear();
  event.theta.clear();
  event.phi.clear();
  event.dp.clear();
  event.path.clear();
  event.m2.clear();
  event.x0s.clear();
  event.y0s.clear();
  event.u0s.clear();
  event.v0s.clear();
  event.poss.clear();
  event.ress.clear();
}


bool EventScat1Tracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventScat1Tracking;
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

  //Generated Scat
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

  //SFT
  tree->Branch("sftnhits", &event.sftnhits);
  tree->Branch("sftlayer", &event.sftlayer);
  tree->Branch("sftposx", &event.sftposx);
  tree->Branch("sftposy", &event.sftposy);
  tree->Branch("sftdl", &event.sftdl);

  //IT1
  tree->Branch("it1nhits", &event.it1nhits);
  tree->Branch("it1layer", &event.it1layer);
  tree->Branch("it1posx", &event.it1posx);
  tree->Branch("it1posy", &event.it1posy);
  tree->Branch("it1dl", &event.it1dl);

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

  //Local Tracking
  tree->Branch("ntslin1", &event.ntslin1);
  tree->Branch("layerslin1", &event.layerslin1);
  tree->Branch("chisqrslin1", &event.chisqrslin1);
  tree->Branch("x0slin1", &event.x0slin1);
  tree->Branch("y0slin1", &event.y0slin1);
  tree->Branch("u0slin1", &event.u0slin1);
  tree->Branch("v0slin1", &event.v0slin1);
  tree->Branch("posslin1", &event.posslin1);
  tree->Branch("resslin1", &event.resslin1);

  tree->Branch("ntslin2", &event.ntslin2);
  tree->Branch("layerslin2", &event.layerslin2);
  tree->Branch("chisqrslin2", &event.chisqrslin2);
  tree->Branch("x0slin2", &event.x0slin2);
  tree->Branch("y0slin2", &event.y0slin2);
  tree->Branch("u0slin2", &event.u0slin2);
  tree->Branch("v0slin2", &event.v0slin2);
  tree->Branch("posslin2", &event.posslin2);
  tree->Branch("resslin2", &event.resslin2);

  tree->Branch("ntslout1", &event.ntslout1);
  tree->Branch("layerslout1", &event.layerslout1);
  tree->Branch("chisqrslout1", &event.chisqrslout1);
  tree->Branch("x0slout1", &event.x0slout1);
  tree->Branch("y0slout1", &event.y0slout1);
  tree->Branch("u0slout1", &event.u0slout1);
  tree->Branch("v0slout1", &event.v0slout1);
  tree->Branch("posslout1", &event.posslout1);
  tree->Branch("resslout1", &event.resslout1);

  tree->Branch("ntslout2", &event.ntslout2);
  tree->Branch("layerslout2", &event.layerslout2);
  tree->Branch("chisqrslout2", &event.chisqrslout2);
  tree->Branch("x0slout2", &event.x0slout2);
  tree->Branch("y0slout2", &event.y0slout2);
  tree->Branch("u0slout2", &event.u0slout2);
  tree->Branch("v0slout2", &event.v0slout2);
  tree->Branch("posslout2", &event.posslout2);
  tree->Branch("resslout2", &event.resslout2);

  //Scat Tracking
  tree->Branch("ids", &event.ids);
  tree->Branch("types", &event.types);
  tree->Branch("ntsprein", &event.ntsprein);
  tree->Branch("ntspreout", &event.ntspreout);
  tree->Branch("nts", &event.nts);
  tree->Branch("layers", &event.layers);
  tree->Branch("chisqrs", &event.chisqrs);
  tree->Branch("p", &event.p);
  tree->Branch("theta", &event.theta);
  tree->Branch("phi", &event.phi);
  tree->Branch("dp", &event.dp);
  tree->Branch("path", &event.path);
  tree->Branch("m2", &event.m2);
  tree->Branch("x0s", &event.x0s);
  tree->Branch("y0s", &event.y0s);
  tree->Branch("u0s", &event.u0s);
  tree->Branch("v0s", &event.v0s);
  tree->Branch("poss", &event.poss);
  tree->Branch("ress", &event.ress);

  return true;
}
