/*
  UserScat3Tracking.cc

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

class EventScat3Tracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  HodoAnalyzer *HodoAna;

public:
  EventScat3Tracking();
  ~EventScat3Tracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventScat3Tracking::EventScat3Tracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer()),
    HodoAna(new HodoAnalyzer())
{
}

EventScat3Tracking::~EventScat3Tracking()
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
  int    sftnhits;
  std::vector<int> sftlayer;
  std::vector<double> sftposx, sftposy;
  std::vector<double> sftdl;

  //IT1
  int    it1nhits;
  std::vector<int> it1layer;
  std::vector<double> it1posx, it1posy;
  std::vector<double> it1dl;

  //PAD
  int    padnhits;
  std::vector<int>    padlayer, padseg;
  std::vector<double> padtime, padedep;
  std::vector<double> padpath, padp;
  std::vector<double> padposx, padposy;
  std::vector<int>    padpid;
  std::vector<double> padbeta;

  //Local tracking 
  int    ntslin1;
  std::vector<int>    layerslin1;
  std::vector<double> chisqrslin1;
  std::vector<double> x0slin1, y0slin1;
  std::vector<double> u0slin1, v0slin1;
  std::vector<double> posslin1, resslin1;

  int    ntslin2;
  std::vector<int>    layerslin2;
  std::vector<double> chisqrslin2;
  std::vector<double> x0slin2, y0slin2;
  std::vector<double> u0slin2, v0slin2;
  std::vector<double> posslin2, resslin2;

  //Scat Tracking
  int    nts;
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

bool EventScat3Tracking::ProcessingBegin()
{
 return true;
}

bool EventScat3Tracking::ProcessingNormal( std::ifstream &In )
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

      //PAD
      const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
      int nhPAD = contPAD.size();
      event.padnhits = nhPAD;
      for(int j=0; j<nhPAD; j++){
       	s_HodoRawHit *hitPAD=contPAD[j];
	event.padlayer.push_back(hitPAD->LayerId());
      	event.padseg.push_back(hitPAD->SegmentId());

	//std::cout<< "padlayer= " << hitPAD->LayerId() << std::endl;

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

      //SFT
      for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      	const s_TrRHitContainer &contSFT = hit->GetsSFTRHC(layer);
      	int nhSFT=contSFT.size();
      	event.sftnhits = nhSFT;
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
      	event.it1nhits = nhIT1;
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

      //////////////////////////Tracking
      event.ids.push_back(hit->TrackId());
      event.types.push_back(hit->TrackType());
      IniP=hit->GetMom();

      //if( hit->TrackType() == TrackTypeScat3 ){
	//std::cout<< "*****Before Tracking" << std::endl;
	TrAna->DecodesSRawHits( hit, hit->TrackType() );

	// //ScatInT
	// TrAna->TrackSearchScatInT();
	// int ntin1=TrAna->GetNtracksScatInT();
	// event.ntslin1=ntin1;
	// for( int it=0; it<ntin1; ++it ){
	//   TrLocalTrack *tp=TrAna->GetTrackScatInT(it);
	//   int nh=tp->GetNHit();
	//   double chisqr=tp->GetChiSquare();
	//   double x0=tp->GetX0(), y0=tp->GetY0();
	//   double u0=tp->GetU0(), v0=tp->GetV0();
	  
	//   int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
	//   double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
	//   double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	//   double utgt=u0, vtgt=v0;
	  
	//   event.chisqrslin1.push_back(chisqr);
	//   event.x0slin1.push_back(xtgt);
	//   event.y0slin1.push_back(ytgt);
	//   event.u0slin1.push_back(utgt);
	//   event.v0slin1.push_back(vtgt); 
	  
	//   for( int ih=0; ih<nh; ++ih ){
	//     TrLTrackHit *hit=tp->GetHit(ih);
	//     int layerId=hit->GetLayer(); 
	//     event.layerslin1.push_back(layerId);  
	//     double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	//     event.posslin1.push_back(pos);
	//     event.resslin1.push_back(res);
	//   }
	// }

	//SFT
	TrAna->TrackSearchSFTT();
	int ntin1=TrAna->GetNtracksSFTT();
	event.ntslin1=ntin1;
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
	event.ntslin2=ntin2;
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
	
	//Scat3 tracking
      	TrAna->TrackSearchScat3( IniP, IniVert );
      	int ntScat3=TrAna->GetNTracksScat3();
      	event.nts = ntScat3;

      	for( int it=0; it<ntScat3; ++it ){
      	  Scat3Track *tp=TrAna->GetScat3Track(it);
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

	  //PAD
	  const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
	  int nhPAD = contPAD.size();
	  for(int j=0; j<nhPAD; j++){
	    s_HodoRawHit *hitPAD=contPAD[j];
	    
	    int nt = hitPAD->GetSize();
	    for( int k=0; k<nt; k++ ){
	      double time = hitPAD->GetTime(k);
	      double path = hitPAD->GetPath(k);
	      double mom = hitPAD->GetMom(k);
	      
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

	//    }
    }
  }


  tree->Fill();

  return true;
}

void EventScat3Tracking::InitializeEvent( void )
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
  event.sftnhits = -1;
  event.sftlayer.clear();
  event.sftposx.clear();
  event.sftposy.clear();
  event.sftdl.clear();

  //IT1
  event.it1nhits = -1;
  event.it1layer.clear();
  event.it1posx.clear();
  event.it1posy.clear();
  event.it1dl.clear();

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

  //Local Tracking
  event.ntslin1 = -1;
  event.layerslin1.clear();
  event.chisqrslin1.clear();
  event.x0slin1.clear();
  event.y0slin1.clear();
  event.u0slin1.clear();
  event.v0slin1.clear();
  event.posslin1.clear();
  event.resslin1.clear();

  event.ntslin2 = -1;
  event.layerslin2.clear();
  event.chisqrslin2.clear();
  event.x0slin2.clear();
  event.y0slin2.clear();
  event.u0slin2.clear();
  event.v0slin2.clear();
  event.posslin2.clear();
  event.resslin2.clear();

  //Scat Tracking
  event.ids.clear();
  event.types.clear();
  event.nts = -1;
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


bool EventScat3Tracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventScat3Tracking;
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

  //Scat Tracking
  tree->Branch("ids", &event.ids);
  tree->Branch("types", &event.types);
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
