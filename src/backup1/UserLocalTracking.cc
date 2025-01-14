/*
  UserLocalTracking.cc

  2016/3 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "Main.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "PrimInfo.hh"
#include "SpecLib.hh"

#define check 0
#define bssd_tracking 1
#define sssd1_tracking 1
#define sssd2_tracking 1

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventLocalTracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  // HodoAnalyzer *HodoAna;

public:
  EventLocalTracking();
  ~EventLocalTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( TFile* );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventLocalTracking::EventLocalTracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer())
    // HodoAna(new HodoAnalyzer())
{
}

EventLocalTracking::~EventLocalTracking()
{
  // if (HodoAna)   delete HodoAna;
  if (TrAna)   delete TrAna;
  if (rawData) delete rawData;
}

struct Event{
  //T0
  int    t0nhits;
  std::vector<int>    t0segment;
  std::vector<double> t0tdctop, t0tdcbot;
  std::vector<double> t0adctop, t0adcbot;

  //bSSD
  int    bssdnhits;
  std::vector<int>    bssdlayer;
  std::vector<double> bssdposx, bssdposy;
  std::vector<double> bssddl;

  //bSSD Tracking
  int    ntbssd;
  std::vector<int>    layerbssd;
  std::vector<double> chisqrbssd;
  std::vector<double> x0bssd, y0bssd;
  std::vector<double> u0bssd, v0bssd;
  std::vector<double> posbssd, resbssd;

  //RPC
  int    rpcnhits;
  std::vector<int>    rpcsegment;
  std::vector<double> rpctdclft, rpctdcrgt;
  std::vector<double> rpctotlft, rpctotrgt;

  //sSSD1
  int    sssd1nhits;
  std::vector<int>    sssd1layer;
  std::vector<double> sssd1posx, sssd1posy;
  std::vector<double> sssd1dl;

  //sSSD1 Tracking
  int    ntsssd1;
  std::vector<int>    layersssd1;
  std::vector<double> chisqrsssd1;
  std::vector<double> x0sssd1, y0sssd1;
  std::vector<double> u0sssd1, v0sssd1;
  std::vector<double> possssd1, ressssd1;

  //sSSD2
  int    sssd2nhits;
  std::vector<int>    sssd2layer;
  std::vector<double> sssd2posx, sssd2posy;
  std::vector<double> sssd2dl;

  //sSSD2 Tracking
  int    ntsssd2;
  std::vector<int>    layersssd2;
  std::vector<double> chisqrsssd2;
  std::vector<double> x0sssd2, y0sssd2;
  std::vector<double> u0sssd2, v0sssd2;
  std::vector<double> possssd2, ressssd2;

  // //Primary
  // int    prinhits;
  // double priposx, priposy, priposz;
  // double pbeam, ubeam, vbeam, abmom;
  // double prim1, prip1;
  // double pritheta1, priphi1, prithetacm1, priphicm1;
  // double prim2, prip2;
  // double pritheta2, priphi2, prithetacm2, priphicm2;

  // //Generated Scat
  // int    gsnhits;
  // std::vector<int>    gsid, gstype;
  // std::vector<double> gsp, gspx, gspy, gspz;
  // std::vector<int>    gspid;
  // std::vector<double> gsvx, gsvy;

  // //SFT
  // int    sftnhits;
  // std::vector<int> sftlayer;
  // std::vector<double> sftposx, sftposy;
  // std::vector<double> sftdl;

  // //IT1
  // int    it1nhits;
  // std::vector<int> it1layer;
  // std::vector<double> it1posx, it1posy;
  // std::vector<double> it1dl;

  // //PAD
  // int    padnhits;
  // std::vector<int>    padlayer, padseg;
  // std::vector<double> padtime, padedep;
  // std::vector<double> padpath, padp;
  // std::vector<double> padposx, padposy;
  // std::vector<int>    padpid;
  // std::vector<double> padbeta;

  // //Local tracking 
  // int    ntslin1;
  // std::vector<int>    layerslin1;
  // std::vector<double> chisqrslin1;
  // std::vector<double> x0slin1, y0slin1;
  // std::vector<double> u0slin1, v0slin1;
  // std::vector<double> posslin1, resslin1;

  // int    ntslin2;
  // std::vector<int>    layerslin2;
  // std::vector<double> chisqrslin2;
  // std::vector<double> x0slin2, y0slin2;
  // std::vector<double> u0slin2, v0slin2;
  // std::vector<double> posslin2, resslin2;

  // //Scat Tracking
  // int    nts;
  // std::vector<int>    ids, types;
  // std::vector<int>    layers;
  // std::vector<double> chisqrs;
  // std::vector<double> p, theta, phi;
  // std::vector<double> dp;
  // std::vector<double> path, m2;
  // std::vector<double> x0s, y0s;
  // std::vector<double> u0s, v0s;
  // std::vector<double> poss, ress;
};
static Event event;

bool EventLocalTracking::ProcessingBegin()
{
 return true;
}

bool EventLocalTracking::ProcessingNormal( TFile* iFile )
{
  const std::string funcname = "ProcessingNormal";

  TDirectoryFile *iDir = (TDirectoryFile*)iFile->Get("out1");
  TTree *tree = (TTree*)iDir->Get("T0AnaTree");

  int n_entry = tree->GetEntries();
  for(int i_entry = 0; i_entry < n_entry; i_entry++){

    rawData = new RawData;
    if( !rawData->DecodeRawHits(iFile, i_entry) ) return false;

#if check
    std::cout << "*************************************************" << std::endl; 
#endif

    const TrGeomMan & geomMan=TrGeomMan::GetInstance();

    //**************************************************************************
    //******************RawData

    TTree *tree = dynamic_cast<TTree *>(oFile->Get("tree"));
    InitializeEvent();


    const s_BeamRHitContainer &bcont =rawData->GetsBeamRHC();
    int nhB = bcont.size();
    for(int i=0; i<nhB; i++){
      s_BeamRawHit *hit=bcont[i];
      TrAna->DecodesBRawHits( hit, hit->TrackType() );

      //T0
      const s_HodoRHitContainer &contT0 = hit->GetsT0RHC();
      int nhT0=contT0.size();
      event.t0nhits = nhT0;
      for( int j=0; j<nhT0; ++j){
	s_HodoRawHit *hitT0=contT0[j];
	event.t0segment.push_back(hitT0->SegmentId());

	int nt = hitT0->GetSize();
	for( int k=0; k<nt; k++ ) {
	  event.t0tdctop.push_back(hitT0->GetTDC0Time(k));
	  event.t0tdcbot.push_back(hitT0->GetTDC1Time(k));
	  event.t0adctop.push_back(hitT0->GetADC0Hgt(k));
	  event.t0adcbot.push_back(hitT0->GetADC1Hgt(k));
	}//for(k:nt)
      }//for(j:nht0)


      //bSSD
      for( int layer=1; layer<=NumOfLayersbSSD; ++layer ){
	const s_TrRHitContainer &contbSSD = hit->GetsbSSDRHC(layer);
	int nhbSSD=contbSSD.size();
	event.bssdnhits = nhbSSD;
	for( int j=0; j<nhbSSD; ++j ){
	  s_TrRawHit *hitbSSD=contbSSD[j];
	  event.bssdlayer.push_back(hitbSSD->LayerId());

	  int nt = hitbSSD->GetSize();
	  for( int k=0; k<nt; k++ ) {
	    event.bssdposx.push_back(hitbSSD->GetPosX(k));
	    event.bssdposy.push_back(hitbSSD->GetPosY(k));
	    event.bssddl.push_back(hitbSSD->GetDL(k));
	  }//for(k:nt)
	}//for(j:nhbSSD)
      }//for(layer:NumOfLayersbSSD)


#if bssd_tracking
      //bSSD Tracking
      TrAna->TrackSearchbSSDT();
      int ntbssd = TrAna->GetNtracksbSSDT();
      if(ntbssd==-1){
	std::cout << "Some Wrong things happens" << std::endl;
      }
      event.ntbssd = ntbssd;
      for(int it=0; it<ntbssd; it++){
	TrLocalTrack *tp=TrAna->GetTrackbSSDT(it);
	int nh=tp->GetNHit();
	double chisqr=tp->GetChiSquare();
	double x0=tp->GetX0(), y0=tp->GetY0();
	double u0=tp->GetU0(), v0=tp->GetV0();

	int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
	double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
  
	double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	double utgt=u0, vtgt=v0;

	event.chisqrbssd.push_back(chisqr);
	event.x0bssd.push_back(xtgt);
	event.y0bssd.push_back(ytgt);
	event.u0bssd.push_back(utgt);
	event.v0bssd.push_back(vtgt); 

	for( int ih=0; ih<nh; ++ih ){
	  TrLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer(); 
	  event.layerbssd.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.posbssd.push_back(pos);
	  event.resbssd.push_back(res);
	}//for(ih:nh)
      }//for(it:ntbssd)
#endif

    }//for(i:nhB)


    const s_ScatRHitContainer &scont =rawData->GetsScatRHC();
    int nhS = scont.size();
    for(int i=0; i<nhS; i++){
      s_ScatRawHit *hit=scont[i];
      TrAna->DecodesSRawHits( hit, hit->TrackType() );

      //RPC
      const s_HodoRHitContainer &contRPC = hit->GetsRPCRHC();
      int nhRPC=contRPC.size();
      event.rpcnhits = nhRPC;
      for( int j=0; j<nhRPC; ++j){
	s_HodoRawHit *hitRPC=contRPC[j];
	event.rpcsegment.push_back(hitRPC->SegmentId());

	int nt = hitRPC->GetSize();
	for( int k=0; k<nt; k++ ) {
	  event.rpctdclft.push_back(hitRPC->GetTDC0Time(k));
	  event.rpctdcrgt.push_back(hitRPC->GetTDC1Time(k));
	  event.rpctotlft.push_back(hitRPC->GetTDC0Tot(k));
	  event.rpctotrgt.push_back(hitRPC->GetTDC1Tot(k));
	}//for(k:nt)
      }//for(j:nhrpc)

      //sSSD1
      for( int layer=1; layer<=NumOfLayerssSSD1; ++layer ){
	const s_TrRHitContainer &contsSSD1 = hit->GetssSSD1RHC(layer);
	int nhsSSD1=contsSSD1.size();
	event.sssd1nhits = nhsSSD1;
	for( int j=0; j<nhsSSD1; ++j ){
	  s_TrRawHit *hitsSSD1=contsSSD1[j];
	  event.sssd1layer.push_back(hitsSSD1->LayerId());

	  int nt = hitsSSD1->GetSize();
	  for( int k=0; k<nt; k++ ) {
	    event.sssd1posx.push_back(hitsSSD1->GetPosX(k));
	    event.sssd1posy.push_back(hitsSSD1->GetPosY(k));
	    event.sssd1dl.push_back(hitsSSD1->GetDL(k));
	  }//for(k:nt)
	}//for(j_nhsSSD1)
      }//for(layer:NumOfLayerssSSD1)


#if sssd1_tracking
      //sSSD1 Tracking
      TrAna->TrackSearchsSSD1T();
      int ntsssd1 = TrAna->GetNtrackssSSD1T();
      event.ntsssd1 = ntsssd1;
      for(int it=0; it<ntsssd1; it++){
	TrLocalTrack *tp=TrAna->GetTracksSSD1T(it);
	int nh=tp->GetNHit();
	double chisqr=tp->GetChiSquare();
	double x0=tp->GetX0(), y0=tp->GetY0();
	double u0=tp->GetU0(), v0=tp->GetV0();

	int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
	double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
  
	double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	double utgt=u0, vtgt=v0;

	event.chisqrsssd1.push_back(chisqr);
	event.x0sssd1.push_back(xtgt);
	event.y0sssd1.push_back(ytgt);
	event.u0sssd1.push_back(utgt);
	event.v0sssd1.push_back(vtgt); 

	for( int ih=0; ih<nh; ++ih ){
	  TrLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer(); 
	  event.layersssd1.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.possssd1.push_back(pos);
	  event.ressssd1.push_back(res);
	}//for(ih:nh)
      }//for(it:ntsssd1)
#endif


      //sSSD2
      for( int layer=1; layer<=NumOfLayerssSSD2; ++layer ){
	const s_TrRHitContainer &contsSSD2 = hit->GetssSSD2RHC(layer);
	int nhsSSD2=contsSSD2.size();
	event.sssd2nhits = nhsSSD2;
	for( int j=0; j<nhsSSD2; ++j ){
	  s_TrRawHit *hitsSSD2=contsSSD2[j];
	  event.sssd2layer.push_back(hitsSSD2->LayerId());
      
	  int nt = hitsSSD2->GetSize();
	  for( int k=0; k<nt; k++ ) {
	    event.sssd2posx.push_back(hitsSSD2->GetPosX(k));
	    event.sssd2posy.push_back(hitsSSD2->GetPosY(k));
	    event.sssd2dl.push_back(hitsSSD2->GetDL(k));
	  }//for(k:nt)
	}//for(j:nhsSSD2)
      }//for(layer:NumOfLayerssSSD2)


#if sssd2_tracking
      //sSSD2 Tracking
      TrAna->TrackSearchsSSD2T();
      int ntsssd2 = TrAna->GetNtrackssSSD2T();
      event.ntsssd2 = ntsssd2;
      for(int it=0; it<ntsssd2; it++){
	TrLocalTrack *tp=TrAna->GetTracksSSD2T(it);
	int nh=tp->GetNHit();
	double chisqr=tp->GetChiSquare();
	double x0=tp->GetX0(), y0=tp->GetY0();
	double u0=tp->GetU0(), v0=tp->GetV0();

	int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
	double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
  
	double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
	double utgt=u0, vtgt=v0;

	event.chisqrsssd2.push_back(chisqr);
	event.x0sssd2.push_back(xtgt);
	event.y0sssd2.push_back(ytgt);
	event.u0sssd2.push_back(utgt);
	event.v0sssd2.push_back(vtgt); 

	for( int ih=0; ih<nh; ++ih ){
	  TrLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer(); 
	  event.layersssd2.push_back(layerId);  
	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  event.possssd2.push_back(pos);
	  event.ressssd2.push_back(res);
	}//for(ih:nh)
      }//for(it:ntsssd2)
#endif


    }//for(i:nhS)


  // double X0, Y0, Z0;
  // //PrimaryInfo
  // {
  //   const PrimInfoContainer &cont=rawData->GetPrimHC();  
  //   int nh=cont.size();
  //   for( int i=0; i<nh; ++i ){
  //     PrimInfo *hit=cont[i];
  //     event.priposx = hit->GetVertX();
  //     event.priposy = hit->GetVertY();
  //     event.priposz = hit->GetVertZ();
  //     X0=hit->GetVertX();
  //     Y0=hit->GetVertY();
  //     Z0=hit->GetVertZ();
      
  //     event.pbeam = hit->GetBeamMom();
  //     event.ubeam = hit->GetBeamU();
  //     event.vbeam = hit->GetBeamV();
  //     event.abmom = hit->GetAnaBeamMom();
      
  //     event.prim1 = hit->GetMass1();
  //     event.prip1 = hit->GetMom1();
  //     event.pritheta1 = hit->GetTheta1();
  //     event.priphi1 = hit->GetPhi1();
  //     event.prithetacm1 = hit->GetThetaCM1();
  //     event.priphicm1 = hit->GetPhiCM1();
      
  //     event.prim2 = hit->GetMass2();
  //     event.prip2 = hit->GetMom2();
  //     event.pritheta2 = hit->GetTheta2();
  //     event.priphi2 = hit->GetPhi2();
  //     event.prithetacm2 = hit->GetThetaCM2();
  //     event.priphicm2 = hit->GetPhiCM2();
  //   }
  // }  
  // ThreeVector IniVert(1400.0-Z0, 0.0, 0.0);

//   //Scat
//   double IniP;
//   {
//     const s_ScatRHitContainer &cont =rawData->GetsScatRHC();
//     int nhS = cont.size();
//     event.gsnhits = nhS;
//     for(int i=0; i<nhS; i++){
//       s_ScatRawHit *hit=cont[i];

//       event.gsid.push_back(hit->TrackId());
//       event.gstype.push_back(hit->TrackType());
//       event.gsp.push_back(hit->GetMom());
//       event.gspx.push_back(hit->GetMomX());
//       event.gspy.push_back(hit->GetMomY());
//       event.gspz.push_back(hit->GetMomZ());
//       event.gspid.push_back(hit->GetPid());
//       event.gsvx.push_back(hit->GetVertX());
//       event.gsvy.push_back(hit->GetVertY());

//       //PAD
//       const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
//       int nhPAD = contPAD.size();
//       event.padnhits = nhPAD;
//       for(int j=0; j<nhPAD; j++){
//        	s_HodoRawHit *hitPAD=contPAD[j];
// 	event.padlayer.push_back(hitPAD->LayerId());
//       	event.padseg.push_back(hitPAD->SegmentId());

// 	//std::cout<< "padlayer= " << hitPAD->LayerId() << std::endl;

//       	int nt = hitPAD->GetSize();
// 	for( int k=0; k<nt; k++ ){
//       	  event.padtime.push_back(hitPAD->GetTime(k));
//       	  event.padedep.push_back(hitPAD->GetEdep(k));
//       	  event.padpath.push_back(hitPAD->GetPath(k));
//       	  event.padp.push_back(hitPAD->GetMom(k));
//       	  event.padposx.push_back(hitPAD->GetPosX(k));
//       	  event.padposy.push_back(hitPAD->GetPosY(k));
//       	  event.padpid.push_back(hitPAD->GetPid(k));
//       	  event.padbeta.push_back(hitPAD->GetBeta(k));
// 	}
//       }

//       //SFT
//       for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
//       	const s_TrRHitContainer &contSFT = hit->GetsSFTRHC(layer);
//       	int nhSFT=contSFT.size();
//       	event.sftnhits = nhSFT;
//       	for( int j=0; j<nhSFT; ++j ){
//       	  s_TrRawHit *hitSFT=contSFT[j];
//       	  event.sftlayer.push_back(hitSFT->LayerId());

//       	  int nt = hitSFT->GetSize();
//       	  for( int k=0; k<nt; k++ ) {
//       	    event.sftposx.push_back(hitSFT->GetPosX(k));
//       	    event.sftposy.push_back(hitSFT->GetPosY(k));
//       	    event.sftdl.push_back(hitSFT->GetDL(k));
	
//       	  }
//       	}
//       }

//       //IT1
//       for( int layer=1; layer<=NumOfLayersIT1; ++layer ){
//       	const s_TrRHitContainer &contIT1 = hit->GetsIT1RHC(layer);
//       	int nhIT1=contIT1.size();
//       	event.it1nhits = nhIT1;
//       	for( int j=0; j<nhIT1; ++j ){
//       	  s_TrRawHit *hitIT1=contIT1[j];
//       	  event.it1layer.push_back(hitIT1->LayerId());

//       	  int nt = hitIT1->GetSize();
//       	  for( int k=0; k<nt; k++ ) {
//       	    event.it1posx.push_back(hitIT1->GetPosX(k));
//       	    event.it1posy.push_back(hitIT1->GetPosY(k));
//       	    event.it1dl.push_back(hitIT1->GetDL(k));
//       	  }
//       	}
//       }

      // //////////////////////////Tracking
      // event.ids.push_back(hit->TrackId());
      // event.types.push_back(hit->TrackType());
      // IniP=hit->GetMom();

//       //if( hit->TrackType() == TrackTypeLocal ){
// 	//std::cout<< "*****Before Tracking" << std::endl;
// 	TrAna->DecodesSRawHits( hit, hit->TrackType() );

// 	// //ScatInT
// 	// TrAna->TrackSearchScatInT();
// 	// int ntin1=TrAna->GetNtracksScatInT();
// 	// event.ntslin1=ntin1;
// 	// for( int it=0; it<ntin1; ++it ){
// 	//   TrLocalTrack *tp=TrAna->GetTrackScatInT(it);
// 	//   int nh=tp->GetNHit();
// 	//   double chisqr=tp->GetChiSquare();
// 	//   double x0=tp->GetX0(), y0=tp->GetY0();
// 	//   double u0=tp->GetU0(), v0=tp->GetV0();
	  
// 	//   int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
// 	//   double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
// 	//   double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
// 	//   double utgt=u0, vtgt=v0;
	  
// 	//   event.chisqrslin1.push_back(chisqr);
// 	//   event.x0slin1.push_back(xtgt);
// 	//   event.y0slin1.push_back(ytgt);
// 	//   event.u0slin1.push_back(utgt);
// 	//   event.v0slin1.push_back(vtgt); 
	  
// 	//   for( int ih=0; ih<nh; ++ih ){
// 	//     TrLTrackHit *hit=tp->GetHit(ih);
// 	//     int layerId=hit->GetLayer(); 
// 	//     event.layerslin1.push_back(layerId);  
// 	//     double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
// 	//     event.posslin1.push_back(pos);
// 	//     event.resslin1.push_back(res);
// 	//   }
// 	// }

// 	//SFT
// 	TrAna->TrackSearchSFTT();
// 	int ntin1=TrAna->GetNtracksSFTT();
// 	event.ntslin1=ntin1;
// 	for( int it=0; it<ntin1; ++it ){
// 	  TrLocalTrack *tp=TrAna->GetTrackSFTT(it);
// 	  int nh=tp->GetNHit();
// 	  double chisqr=tp->GetChiSquare();
// 	  double x0=tp->GetX0(), y0=tp->GetY0();
// 	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
// 	  int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
// 	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
// 	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
// 	  double utgt=u0, vtgt=v0;
	  
// 	  event.chisqrslin1.push_back(chisqr);
// 	  event.x0slin1.push_back(xtgt);
// 	  event.y0slin1.push_back(ytgt);
// 	  event.u0slin1.push_back(utgt);
// 	  event.v0slin1.push_back(vtgt); 
	  
// 	  for( int ih=0; ih<nh; ++ih ){
// 	    TrLTrackHit *hit=tp->GetHit(ih);
// 	    int layerId=hit->GetLayer(); 
// 	    event.layerslin1.push_back(layerId);  
// 	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
// 	    event.posslin1.push_back(pos);
// 	    event.resslin1.push_back(res);
// 	  }
// 	}
	
// 	//IT1
// 	TrAna->TrackSearchIT1T();
// 	int ntin2=TrAna->GetNtracksIT1T();
// 	event.ntslin2=ntin2;
// 	for( int it=0; it<ntin2; ++it ){
// 	  TrLocalTrack *tp=TrAna->GetTrackIT1T(it);
// 	  int nh=tp->GetNHit();
// 	  double chisqr=tp->GetChiSquare();
// 	  double x0=tp->GetX0(), y0=tp->GetY0();
// 	  double u0=tp->GetU0(), v0=tp->GetV0();
	  
// 	  int Id = TrGeomMan::GetInstance().GetDetectorId("Target");
// 	  double z = TrGeomMan::GetInstance().GetLocalZ( Id ); 
	  
// 	  double xtgt=tp->GetX( z ), ytgt=tp->GetY( z );
// 	  double utgt=u0, vtgt=v0;
	  
// 	  event.chisqrslin2.push_back(chisqr);
// 	  event.x0slin2.push_back(xtgt);
// 	  event.y0slin2.push_back(ytgt);
// 	  event.u0slin2.push_back(utgt);
// 	  event.v0slin2.push_back(vtgt); 
	  
// 	  for( int ih=0; ih<nh; ++ih ){
// 	    TrLTrackHit *hit=tp->GetHit(ih);
// 	    int layerId=hit->GetLayer(); 
// 	    event.layerslin2.push_back(layerId);  
// 	    double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
// 	    event.posslin2.push_back(pos);
// 	    event.resslin2.push_back(res);
// 	  }
// 	}
	
// //       	TrAna->TrackSearchScat3( IniP, IniVert );
// //       	int ntScat3=TrAna->GetNTracksScat3();
// //       	event.nts = ntScat3;

// //       	for( int it=0; it<ntScat3; ++it ){
// //       	  Scat3Track *tp=TrAna->GetScat3Track(it);
// //       	  if(!tp) continue;
// //       	  int nh=tp->GetNHits();
// //       	  double chisqr=tp->chisqr();  
// //       	  ThreeVector Ppos=tp->PrimaryPosition();
// //       	  ThreeVector Pmom=tp->PrimaryMomentum();
// //       	  double pathL=tp->PathLengthToTOF();
// //       	  double x=Ppos.x(), y=Ppos.y();
// //       	  double p=Pmom.mag();
// // 	  double pz=Pmom.z();
// //       	  double u, v;
// //       	  double theta, phi;	  

// // 	  if(pz>0){
// // 	    u=Pmom.x()/pz, v=Pmom.y()/pz;
// // 	    theta = Pmom.theta();	  
// // 	    phi = Pmom.phi();	  
// // 	  }
// // 	  if(pz<0){
// // 	    u=Pmom.x()/pz, v=Pmom.y()/pz;
// // 	    theta = (-1.*Pmom).theta();	  
// // 	    phi = (-1.*Pmom).phi();	  
// // 	  }

// //       	  // std::cout<< "V1= " << ThreeVector(X0, Y0, 0.0) << std::endl;
// //       	  // std::cout<< "V2= " << Ppos << std::endl;
	  
// //       	  event.chisqrs.push_back(chisqr);
// //       	  event.p.push_back(p);
// // 	  event.theta.push_back(theta*Rad2Deg);
// // 	  event.phi.push_back(phi*Rad2Deg);
// //       	  event.path.push_back(pathL);
// //       	  event.x0s.push_back(x);
// //       	  event.y0s.push_back(y);
// //       	  event.u0s.push_back(u);
// //       	  event.v0s.push_back(v);   

// // 	  event.dp.push_back(p-IniP);

// // #if check
// // 	  std::cout<< "***** p= " << p << ", dp= " << p-IniP <<std::endl;
// // #endif

// // 	  //PAD
// // 	  const s_HodoRHitContainer &contPAD =hit->GetsPADRHC();
// // 	  int nhPAD = contPAD.size();
// // 	  for(int j=0; j<nhPAD; j++){
// // 	    s_HodoRawHit *hitPAD=contPAD[j];
	    
// // 	    int nt = hitPAD->GetSize();
// // 	    for( int k=0; k<nt; k++ ){
// // 	      double time = hitPAD->GetTime(k);
// // 	      double path = hitPAD->GetPath(k);
// // 	      double mom = hitPAD->GetMom(k);
	      
// // 	      double C = 2.99792458E+8;
// // 	      double T0, V0, B0, m2_0;
	      
// // 	      V0 = (path/1000.)/(time*1.0E-9);
// // 	      B0 = V0/C;
// // 	      m2_0 = (mom/B0)*(mom/B0)*(1.-B0*B0);
	      
// // 	      event.m2.push_back(m2_0); 
// // 	    }
// // 	  }
	  
// // 	  for( int j=0; j<nh; ++j ){
// // 	    TrackHit *hit=tp->GetHit(j);
// // 	    if(!hit) continue;
// // 	    int layerId=hit->GetLayer();
// //       	    event.layers.push_back(layerId);  
// //       	  }
// // 	}

// 	//    }
  //   }
  // }


    tree->Fill();
  }//for(i_entry:n_entry)

  return true;
}

void EventLocalTracking::InitializeEvent( void )
{
  //T0
  event.t0nhits = -1;
  event.t0segment.clear();
  event.t0tdctop.clear();
  event.t0tdcbot.clear();
  event.t0adctop.clear();
  event.t0adcbot.clear();

  //bSSD
  event.bssdnhits = -1;
  event.bssdlayer.clear();
  event.bssdposx.clear();
  event.bssdposy.clear();
  event.bssddl.clear();

  //bSSD Tracking
  event.ntbssd = -1;
  event.layerbssd.clear();
  event.chisqrbssd.clear();
  event.x0bssd.clear();
  event.y0bssd.clear();
  event.u0bssd.clear();
  event.v0bssd.clear();
  event.posbssd.clear();
  event.resbssd.clear();

  //RPC
  event.rpcnhits = -1;
  event.rpcsegment.clear();
  event.rpctdclft.clear();
  event.rpctdcrgt.clear();
  event.rpctotlft.clear();
  event.rpctotrgt.clear();

  //sSSD1
  event.sssd1nhits = -1;
  event.sssd1layer.clear();
  event.sssd1posx.clear();
  event.sssd1posy.clear();
  event.sssd1dl.clear();

  //sSSD1 Tracking
  event.ntsssd1 = -1;
  event.layersssd1.clear();
  event.chisqrsssd1.clear();
  event.x0sssd1.clear();
  event.y0sssd1.clear();
  event.u0sssd1.clear();
  event.v0sssd1.clear();
  event.possssd1.clear();
  event.ressssd1.clear();

  //sSSD2
  event.sssd2nhits = -1;
  event.sssd2layer.clear();
  event.sssd2posx.clear();
  event.sssd2posy.clear();
  event.sssd2dl.clear();

  //sSSD2 Tracking
  event.ntsssd2 = -1;
  event.layersssd2.clear();
  event.chisqrsssd2.clear();
  event.x0sssd2.clear();
  event.y0sssd2.clear();
  event.u0sssd2.clear();
  event.v0sssd2.clear();
  event.possssd2.clear();
  event.ressssd2.clear();

  // //PriInfo
  // event.priposx = -9999.0;
  // event.priposy = -9999.0;
  // event.priposz = -9999.0;
  // event.pbeam = -9999.0;
  // event.ubeam = -9999.0;
  // event.vbeam = -9999.0;
  // event.abmom = -9999.0;
  // event.prim1 = -9999.0;
  // event.prip1 = -9999.0;
  // event.pritheta1 = -9999.0;
  // event.priphi1 = -9999.0;
  // event.prithetacm1 = -9999.0;
  // event.priphicm1 = -9999.0;
  // event.prim2 = -9999.0;
  // event.prip2 = -9999.0;
  // event.pritheta2 = -9999.0;
  // event.priphi2 = -9999.0;
  // event.prithetacm2 = -9999.0;
  // event.priphicm2 = -9999.0;

  // //Generated Scat
  // event.gsnhits = -1;
  // event.gsid.clear();
  // event.gstype.clear();
  // event.gsp.clear();
  // event.gspx.clear();
  // event.gspy.clear();
  // event.gspz.clear();
  // event.gspid.clear();
  // event.gsvx.clear();
  // event.gsvy.clear();

  // //SFT
  // event.sftnhits = -1;
  // event.sftlayer.clear();
  // event.sftposx.clear();
  // event.sftposy.clear();
  // event.sftdl.clear();

  // //IT1
  // event.it1nhits = -1;
  // event.it1layer.clear();
  // event.it1posx.clear();
  // event.it1posy.clear();
  // event.it1dl.clear();

  // //PAD
  // event.padnhits = -1;
  // event.padlayer.clear();
  // event.padseg.clear();
  // event.padtime.clear();
  // event.padedep.clear();
  // event.padpath.clear();
  // event.padp.clear();
  // event.padposx.clear();
  // event.padposy.clear();
  // event.padpid.clear();
  // event.padbeta.clear();

  // //Local Tracking
  // event.ntslin1 = -1;
  // event.layerslin1.clear();
  // event.chisqrslin1.clear();
  // event.x0slin1.clear();
  // event.y0slin1.clear();
  // event.u0slin1.clear();
  // event.v0slin1.clear();
  // event.posslin1.clear();
  // event.resslin1.clear();

  // event.ntslin2 = -1;
  // event.layerslin2.clear();
  // event.chisqrslin2.clear();
  // event.x0slin2.clear();
  // event.y0slin2.clear();
  // event.u0slin2.clear();
  // event.v0slin2.clear();
  // event.posslin2.clear();
  // event.resslin2.clear();

  // //Scat Tracking
  // event.ids.clear();
  // event.types.clear();
  // event.nts = -1;
  // event.layers.clear();
  // event.chisqrs.clear();
  // event.p.clear();
  // event.theta.clear();
  // event.phi.clear();
  // event.dp.clear();
  // event.path.clear();
  // event.m2.clear();
  // event.x0s.clear();
  // event.y0s.clear();
  // event.u0s.clear();
  // event.v0s.clear();
  // event.poss.clear();
  // event.ress.clear();
}


bool EventLocalTracking::ProcessingEnd()
{
  // oFile->Write();
  // oFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventLocalTracking;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(oFile->Get("tree"));

  //T0
  tree->Branch("t0nhits", &event.t0nhits);
  tree->Branch("t0segment", &event.t0segment);
  tree->Branch("t0tdctop", &event.t0tdctop);
  tree->Branch("t0tdcbot", &event.t0tdcbot);
  tree->Branch("t0adctop", &event.t0adctop);
  tree->Branch("t0adcbot", &event.t0adcbot);

  //bSSD
  tree->Branch("bssdnhits", &event.bssdnhits);
  tree->Branch("bssdlayer", &event.bssdlayer);
  tree->Branch("bssdposx", &event.bssdposx);
  tree->Branch("bssdposy", &event.bssdposy);
  tree->Branch("bssddl", &event.bssddl);

  //bSSD Tracking
  tree->Branch("ntbssd", &event.ntbssd);
  tree->Branch("layerbssd", &event.layerbssd);
  tree->Branch("chisqrbssd", &event.chisqrbssd);
  tree->Branch("x0bssd", &event.x0bssd);
  tree->Branch("y0bssd", &event.y0bssd);
  tree->Branch("u0bssd", &event.u0bssd);
  tree->Branch("v0bssd", &event.v0bssd);
  tree->Branch("posbssd", &event.posbssd);
  tree->Branch("resbssd", &event.resbssd);

  //RPC
  tree->Branch("rpcnhits", &event.rpcnhits);
  tree->Branch("rpcsegment", &event.rpcsegment);
  tree->Branch("rpctdclft", &event.rpctdclft);
  tree->Branch("rpctdcrgt", &event.rpctdcrgt);
  tree->Branch("rpctotlft", &event.rpctotlft);
  tree->Branch("rpctotrgt", &event.rpctotrgt);

  //sSSD1
  tree->Branch("sssd1nhits", &event.sssd1nhits);
  tree->Branch("sssd1layer", &event.sssd1layer);
  tree->Branch("sssd1posx", &event.sssd1posx);
  tree->Branch("sssd1posy", &event.sssd1posy);
  tree->Branch("sssd1dl", &event.sssd1dl);

  //sSSD1 Tracking
  tree->Branch("ntsssd1", &event.ntsssd1);
  tree->Branch("layersssd1", &event.layersssd1);
  tree->Branch("chisqrsssd1", &event.chisqrsssd1);
  tree->Branch("x0sssd1", &event.x0sssd1);
  tree->Branch("y0sssd1", &event.y0sssd1);
  tree->Branch("u0sssd1", &event.u0sssd1);
  tree->Branch("v0sssd1", &event.v0sssd1);
  tree->Branch("possssd1", &event.possssd1);
  tree->Branch("ressssd1", &event.ressssd1);

  //sSSD2
  tree->Branch("sssd2nhits", &event.sssd2nhits);
  tree->Branch("sssd2layer", &event.sssd2layer);
  tree->Branch("sssd2posx", &event.sssd2posx);
  tree->Branch("sssd2posy", &event.sssd2posy);
  tree->Branch("sssd2dl", &event.sssd2dl);

  //sSSD2 Tracking
  tree->Branch("ntsssd2", &event.ntsssd2);
  tree->Branch("layersssd2", &event.layersssd2);
  tree->Branch("chisqrsssd2", &event.chisqrsssd2);
  tree->Branch("x0sssd2", &event.x0sssd2);
  tree->Branch("y0sssd2", &event.y0sssd2);
  tree->Branch("u0sssd2", &event.u0sssd2);
  tree->Branch("v0sssd2", &event.v0sssd2);
  tree->Branch("possssd2", &event.possssd2);
  tree->Branch("ressssd2", &event.ressssd2);

  // //PriInfo  
  // tree->Branch("priposx", &event.priposx);
  // tree->Branch("priposy", &event.priposy);
  // tree->Branch("priposz", &event.priposz);
  // tree->Branch("pbeam", &event.pbeam);
  // tree->Branch("ubeam", &event.ubeam);
  // tree->Branch("vbeam", &event.vbeam);
  // tree->Branch("abmom", &event.abmom);
  // tree->Branch("prim1", &event.prim1);
  // tree->Branch("prip1", &event.prip1);
  // tree->Branch("pritheta1", &event.pritheta1);
  // tree->Branch("priphi1", &event.priphi1);
  // tree->Branch("prithetacm1", &event.prithetacm1);
  // tree->Branch("priphicm1", &event.priphicm1);
  // tree->Branch("prim2", &event.prim2);
  // tree->Branch("prip2", &event.prip2);
  // tree->Branch("pritheta2", &event.pritheta2);
  // tree->Branch("priphi2", &event.priphi2);
  // tree->Branch("prithetacm2", &event.prithetacm2);
  // tree->Branch("priphicm2", &event.priphicm2);

  // //Generated Scat
  // tree->Branch("gsnhits", &event.gsnhits);
  // tree->Branch("gsid", &event.gsid);
  // tree->Branch("gstype", &event.gstype);
  // tree->Branch("gsp", &event.gsp);
  // tree->Branch("gspx", &event.gspx);
  // tree->Branch("gspy", &event.gspy);
  // tree->Branch("gspz", &event.gspz);
  // tree->Branch("gspid", &event.gspid);
  // tree->Branch("gsvx", &event.gsvx);
  // tree->Branch("gsvy", &event.gsvy);

  // //SFT
  // tree->Branch("sftnhits", &event.sftnhits);
  // tree->Branch("sftlayer", &event.sftlayer);
  // tree->Branch("sftposx", &event.sftposx);
  // tree->Branch("sftposy", &event.sftposy);
  // tree->Branch("sftdl", &event.sftdl);

  // //IT1
  // tree->Branch("it1nhits", &event.it1nhits);
  // tree->Branch("it1layer", &event.it1layer);
  // tree->Branch("it1posx", &event.it1posx);
  // tree->Branch("it1posy", &event.it1posy);
  // tree->Branch("it1dl", &event.it1dl);

  // //PAD
  // tree->Branch("padnhits", &event.padnhits);
  // tree->Branch("padlayer", &event.padlayer);
  // tree->Branch("padseg", &event.padseg);
  // tree->Branch("padtime", &event.padtime);
  // tree->Branch("padedep", &event.padedep);
  // tree->Branch("padpath", &event.padpath);
  // tree->Branch("padp", &event.padp);
  // tree->Branch("padposx", &event.padposx);
  // tree->Branch("padposy", &event.padposy);
  // tree->Branch("padpid", &event.padpid);
  // tree->Branch("padbeta", &event.padbeta);

  // //Local Tracking
  // tree->Branch("ntslin1", &event.ntslin1);
  // tree->Branch("layerslin1", &event.layerslin1);
  // tree->Branch("chisqrslin1", &event.chisqrslin1);
  // tree->Branch("x0slin1", &event.x0slin1);
  // tree->Branch("y0slin1", &event.y0slin1);
  // tree->Branch("u0slin1", &event.u0slin1);
  // tree->Branch("v0slin1", &event.v0slin1);
  // tree->Branch("posslin1", &event.posslin1);
  // tree->Branch("resslin1", &event.resslin1);

  // tree->Branch("ntslin2", &event.ntslin2);
  // tree->Branch("layerslin2", &event.layerslin2);
  // tree->Branch("chisqrslin2", &event.chisqrslin2);
  // tree->Branch("x0slin2", &event.x0slin2);
  // tree->Branch("y0slin2", &event.y0slin2);
  // tree->Branch("u0slin2", &event.u0slin2);
  // tree->Branch("v0slin2", &event.v0slin2);
  // tree->Branch("posslin2", &event.posslin2);
  // tree->Branch("resslin2", &event.resslin2);

  // //Scat Tracking
  // tree->Branch("ids", &event.ids);
  // tree->Branch("types", &event.types);
  // tree->Branch("nts", &event.nts);
  // tree->Branch("layers", &event.layers);
  // tree->Branch("chisqrs", &event.chisqrs);
  // tree->Branch("p", &event.p);
  // tree->Branch("theta", &event.theta);
  // tree->Branch("phi", &event.phi);
  // tree->Branch("dp", &event.dp);
  // tree->Branch("path", &event.path);
  // tree->Branch("m2", &event.m2);
  // tree->Branch("x0s", &event.x0s);
  // tree->Branch("y0s", &event.y0s);
  // tree->Branch("u0s", &event.u0s);
  // tree->Branch("v0s", &event.v0s);
  // tree->Branch("poss", &event.poss);
  // tree->Branch("ress", &event.ress);

  return true;
}
