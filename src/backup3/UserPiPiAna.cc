/*
  UserPiPiAna.cc
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "SksLib.hh"

const int BEAM_TRIG  = 1;
const int PIPI_TRIG  = 6;
const int BH2_TRIG   = 3;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double offset = 15.5;

const double x_off = 10.28;
const double y_off =  3.00;
const double u_off = 0.0012;
const double v_off = 0.0048;

VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventPiPiAna 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventPiPiAna();
  ~EventPiPiAna();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventPiPiAna::EventPiPiAna()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventPiPiAna::~EventPiPiAna()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 20
#endif

//For Tree
struct Event{
  //Trigger
  int trigtype;
  int trigflag[NumOfMisc];

  //Hodoscope
  double bh1HitSeg[MaxHits];
  double bh2HitSeg[MaxHits];
  double tofHitSeg[MaxHits];
  double lcHitSeg[MaxHits];

  double tbh1Hit[MaxHits];
  double tbh2Hit[MaxHits];
  double debh1Hit[MaxHits];
  double debh2Hit[MaxHits];

  double ttofHit[MaxHits];
  double tlcHit[MaxHits];
  double detofHit[MaxHits];
  double delcHit[MaxHits];

  double btof[MaxHits];

  //DC Beam
  int ntBcIn;
  int nhBcIn[MaxHits]; 
  double chisqrBcIn[MaxHits];
  double x0BcIn[MaxHits];
  double y0BcIn[MaxHits];
  double u0BcIn[MaxHits];
  double v0BcIn[MaxHits];
  double xbh1BcIn[MaxHits];
  double ybh1BcIn[MaxHits];

  int ntBcOut;
  int nhBcOut[MaxHits]; 
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits]; 
  double chisqrK18[MaxHits]; 
  double pK18[MaxHits]; 
  double xtb[MaxHits];
  double ytb[MaxHits];
  double utb[MaxHits];
  double vtb[MaxHits];
  double thetab[MaxHits];

  //DC SKS
  int ntSdcIn;
  int nhSdcIn[MaxHits]; 
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nhSdcOut[MaxHits]; 
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double xtofSdcOut[MaxHits];
  double ytofSdcOut[MaxHits];

  int    ntSks;
  int    nhSks[MaxHits]; 
  double chisqrSks[MaxHits]; 
  double path[MaxHits]; 
  double pSks[MaxHits]; 
  double m2[MaxHits]; 
  double xts[MaxHits];
  double yts[MaxHits];
  double uts[MaxHits];
  double vts[MaxHits];
  double thetas[MaxHits];

  //Reaction
  int    nPi;
  int    nK;
  int    nPiK;
  double vtx[MaxHits];
  double vty[MaxHits];
  double vtz[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double BE[MaxHits];
};
static Event event;

const double MinDeltaEBH2 = 0.2; //for Pion
const double MaxDeltaEBH2 = 4.0; //for Pion

const double MinDeltaEBH1 = 0.2; //for Pion
const double MaxDeltaEBH1 = 4.0; //for Pion

const double MinBeamToF  = -3.0;
const double MaxBeamToF  =  3.0;

const double MinDeltaETof = 0.2; //for Proton
const double MaxDeltaETof = 5.0; //for Proton

const double MinTimeTof = -2.0;
const double MaxTimeTof = 40.0;

//Multi cut of DCs
const int MaxMultiHitBcIn   = 3;
const int MaxMultiHitBcOut  = 100;
const int MaxMultiHitSdcIn  = 100;
const int MaxMultiHitSdcOut = 3;

const double MaxChisqrBcIn   = 100.;
const double MaxChisqrBcOut  = 100.;
const double MaxChisqrSdcIn  = 100.;
const double MaxChisqrSdcOut = 100.;
const double MaxChisqrK18    = 100.;
const double MaxChisqrSks    = 500.;

const int MaxNTrackBcIn   = 3;
const int MaxNTrackBcOut  = 3;
const int MaxNTrackSdcIn  = 3;
const int MaxNTrackSdcOut = 2;
const int MaxNTrackK18    = 2;
const int MaxNTrackSks    = 2;

const double MinTgtXSdcIn  =  -100.;
const double MaxTgtXSdcIn  =   100.;
const double MinTgtYSdcIn  =  -50.;
const double MaxTgtYSdcIn  =   50.;

const double MinTgtXBcOut  =  -100.;
const double MaxTgtXBcOut  =   100.;
const double MinTgtYBcOut  =  -50.;
const double MaxTgtYBcOut  =   50.;

const double MinMassSquare = -0.6;
const double MaxMassSquare = 1.20;

const double MinVertexZ = -20.0;
const double MaxVertexZ =  20.0;

const double AtomicMassUnit  = 0.93149432;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double CarbonMass      = 12.*AtomicMassUnit;
const double DeutronMass     = 2.*AtomicMassUnit+0.01313672;
const double DeltaC11        = 0.0106502;
const double CoreMass        = 11.*AtomicMassUnit+DeltaC11;
const double LambdaMass      = 1.115648;

bool EventPiPiAna::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventPiPiAna::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();
  
  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //Misc
  int trig_type=0;
  int trig=0;
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();
      
      //pi beam
      if( seg==1 && T>0 ){
	trig_type =(trig_type | (1 << 0));
	event.trigtype = (trig_type >> 1)+1;
	trig = (trig_type >> 1)+1;
      }//p beam
      else if( seg==2 && T>0 ){
	trig_type =(trig_type | (1 << 1));
	event.trigtype = (trig_type >> 2)+2;
	trig = (trig_type >> 2)+2;
      }//K baem
      else if( seg==3 && T>0 ){
	trig_type =(trig_type | (1 << 2));
	event.trigtype = (trig_type >> 3)+3;
	trig = (trig_type >> 3)+3;
      }//(pi, K)
      else if( seg==5 && T>0 ){
	trig_type =(trig_type | (1 << 4));
	event.trigtype = (trig_type >> 5)+5;
	trig = (trig_type >> 5)+5;
      }//(pi, pi)
      else if( seg==6 && T>0 ){
	trig_type =(trig_type | (1 << 5));
	event.trigtype = (trig_type >> 6)+6;
	trig = (trig_type >> 6)+6;
      }//(pi, p)
      else if( seg==7 && T>0 ){
	trig_type =(trig_type | (1 << 6));
	event.trigtype = (trig_type >> 7)+7;
	trig = (trig_type >> 7)+7;
      }
      event.trigflag[i] = T;
    }
  }

  //------------------------Cut
  //if (trig != PIPI_TRIG)  return true;
  HF1( 1, 1. );

  ////////////////////////////Hodoscope
  //BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  HF1( 101, double(ncBh2) );
  //------------------------Cut
  if( ncBh2==0 ) {
    //     //EvDisp
    //     if (FlagEvDisp) {
    //       std::cout << "ncBH2=0" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }

  BH2Cluster *clBH2Time0=hodoAna->GetClusterBH2(0);
  double time0=clBH2Time0->CTime0();
  {
    double mint=clBH2Time0->CMeanTime();
    HF1( 102, clBH2Time0->ClusterSize() );
    HF1( 103, clBH2Time0->MeanSeg()-0.5 );
    HF1( 104, mint );
    HF1( 105, clBH2Time0->DeltaE() );
    for( int i=1; i<ncBh2; ++i ){
      BH2Cluster *cl=hodoAna->GetClusterBH2(i);
      double t=cl->CMeanTime();
      double dEbh2 = cl->DeltaE();
      if( !(MinDeltaEBH2<dEbh2 && dEbh2<MaxDeltaEBH2) ) continue;

      HF1( 102, cl->ClusterSize() );
      HF1( 103, cl->MeanSeg()-0.5 );
      HF1( 104, t );
      HF1( 105, cl->DeltaE() );
      if( fabs(t)<fabs(mint) ){
	clBH2Time0=cl;
	mint=t; time0=clBH2Time0->CTime0();
      }
    }
  }
  HF1( 112, clBH2Time0->ClusterSize() );
  HF1( 113, clBH2Time0->MeanSeg()-0.5 );
  HF1( 114, clBH2Time0->CMeanTime() );
  HF1( 115, clBH2Time0->DeltaE() );

  HF1( 1, 2. );
  
  //BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  //------------------------Cut
  {
    int ncOk=0;
    HF1( 201, double(ncBh1) );
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      if( !(MinDeltaEBH1<dEbh1 && dEbh1<MaxDeltaEBH1) ) continue;

      HF1( 202, cl->ClusterSize() );
      HF1( 203, cl->MeanSeg()-0.5 );
      HF1( 204, cl->CMeanTime() );
      HF1( 205, cl->DeltaE() );
      HF1( 206, btof );
      event.btof[ncBh1] = btof;
      //------------------------Cut
      if( btof>MinBeamToF && btof<MaxBeamToF ){
	++ncOk;
	HF1( 212, cl->ClusterSize() );
	HF1( 213, cl->MeanSeg()-0.5 );
	HF1( 214, cl->CMeanTime() );
	HF1( 215, cl->DeltaE() );
      }
      else{
	cl->GoodForAnalysis( false );
      }
    }
    HF1( 211, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //       if (FlagEvDisp) {
      // 	std::cout << "rejected by BH1 cut" << std::endl;
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.get_command();
      // 	evDisp.EndOfEvent();
      //       }
      return true;
    }
  }
  HF1( 1, 3. );

  // Event Display
  //   if (FlagEvDisp){
  //     const HodoRHitContainer &cont=rawData->GetTOFRawHC();
  //     int nh=cont.size();
  //     int TofId=DCGeomMan::GetInstance().GetTofId(); 
  //     for( int i=0; i<nh; ++i ){
  //       HodoRawHit *hit=cont[i];
  //       int seg=hit->SegmentId();
  //       int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
  //       if (Tu>0 || Td>0) {
  // 	if (FlagEvDisp)	{
  // 	  const EvDisp & evDisp = EvDisp::GetInstance();	  
  // 	  evDisp.DrawHitHodoscope(TofId, seg, Tu, Td);
  // 	}
  //       }
  //     }
  //   }

 //Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  HF1( 301, double(ncTof) );
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0;
      HF1( 302, cl->ClusterSize() );
      HF1( 303, cl->MeanSeg()-0.5 );
      HF1( 304, t );
      HF1( 305, de );
      //------------------------Cut
      if( de>MinDeltaETof && de<MaxDeltaETof &&
	  t>MinTimeTof && t<MaxTimeTof ){
	++ncOk;
	HF1( 312, cl->ClusterSize() );
	HF1( 313, cl->MeanSeg()-0.5 );
	HF1( 314, t );
	HF1( 315, de );
      }	
      else{
	cl->GoodForAnalysis( false );
      }
    }
    HF1( 311, double(ncOk) );
    //------------------------Cut
    if( ncOk<1 ) {
      //       if (FlagEvDisp) {
      // 	std::cout << "rejected by Tof cut" << std::endl;
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.get_command();
      // 	evDisp.EndOfEvent();
      //       }
      return true;
    }
  }
  HF1( 1, 4. );

  HF1( 1, 10. );
  ////////////////////////////DC
  int multi_BcIn[12];
  for (int i=0; i<NumOfLayersBcIn; i++) multi_BcIn[i]=0;
  int multi_BcOut[12];
  for (int i=0; i<NumOfLayersBcOut; i++) multi_BcOut[i]=0;
  int multi_SdcIn[10];
  for (int i=0; i<NumOfLayersSdcIn; i++) multi_SdcIn[i]=0;
  int multi_SdcOut[12];
  for (int i=0; i<NumOfLayersSdcOut; i++) multi_SdcOut[i]=0;

  DCAna->DecodeRawHits( rawData ); 
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );
  double zK18Tgt = DCGeomMan::GetInstance().GetLocalZ( IdK18Target );

  //////////////BC1&2 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcInHC(layer);
      int nhIn=contIn.size();
      multi_BcIn[layer-1] = nhIn;
    }
  }
  //------------------------Cut
  for (int i=0; i<NumOfLayersBcIn; i++) {
    if( multi_BcIn[i] >MaxMultiHitBcIn )
      return true;
  }
  
  //////////////BC3&4 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut[layer-1] = nhOut;
    }
  }
  //------------------------Cut
  for (int i=0; i<NumOfLayersBcOut; i++) {
    if( multi_BcOut[i] >MaxMultiHitBcOut )
      return true;
  }
  
  //////////////SDC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn[layer-1] = nhIn;
    }
  }
  //------------------------Cut
  for (int i=0; i<NumOfLayersSdcIn; i++) {
    if( multi_SdcIn[i] >MaxMultiHitSdcIn )
      return true;
  }
  
  //////////////SDC3&4 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut[layer-1] = nhOut;
    }
  }
  //------------------------Cut
  for (int i=0; i<NumOfLayersSdcOut; i++) {
    if( multi_SdcOut[i] >MaxMultiHitSdcOut ){
      if( i!=4 ){
	return true;
      }
    }
  }
  HF1( 1, 11. );

  //SdcIn Analysis
  DCAna->TrackSearchSdcIn();
  int ntSdcIn=DCAna->GetNtracksSdcIn();
  event.ntSdcIn=ntSdcIn;
  HF1( 1001, double(ntSdcIn) );
  //------------------------Cut
  if( ntSdcIn<1 || ntSdcIn>MaxNTrackSdcIn ) {
    //     if (FlagEvDisp) {
    //       std::cout << "rejected by SdcIn" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 12. );

  {
    int ntOk=0;
    for( int i=0; i<ntSdcIn; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      
      if(i<MaxHits){
	event.nhSdcIn[i]=nh;
	event.chisqrSdcIn[i]=chisqr;
	event.u0SdcIn[i]=u0;
	event.v0SdcIn[i]=v0;
	event.x0SdcIn[i]=x0;
	event.y0SdcIn[i]=y0;
      }
      
      HF1( 1002, double(nh) );
      HF1( 1003, chisqr );
      HF1( 1004, x0 ); HF1( 1005, y0 );
      HF1( 1006, u0 ); HF1( 1007, v0 );
      HF2( 1008, x0, u0 ); HF2( 1009, y0, v0 );
      HF2( 1010, x0, y0 );

      //evDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.DrawSdcInLocalTrack(tp);
      //       }
      //------------------------Cut
      if( MinTgtXSdcIn<x0 && x0<MaxTgtXSdcIn && 
	  MinTgtYSdcIn<y0 && y0<MaxTgtYSdcIn &&
	  chisqr<MaxChisqrSdcIn ){
	++ntOk;
      }
    }
    if( ntOk<1 ) {
      //     if (FlagEvDisp) {
      //       std::cout << "rejected by SdcIn" << std::endl;
      //       const EvDisp & evDisp = EvDisp::GetInstance();
      //       evDisp.get_command();
      //       evDisp.EndOfEvent();
      //     }
      return true;
    }
  }
  HF1( 1, 13. );

  //BcOut Analysis
  DCAna->TrackSearchBcOut();
  int ntBcOut=DCAna->GetNtracksBcOut();
  event.ntBcOut=ntBcOut;
  HF1( 1101, double(ntBcOut) );
  //------------------------Cut
  if( ntBcOut<1 || ntBcOut>MaxNTrackBcOut ) {
    //     if (FlagEvDisp) {
    //       std::cout << "rejected by BcOut" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 14. );

  {
    int ntOk=0;
    for( int i=0; i<ntBcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackBcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      double xtgt=tp->GetX( zK18Tgt ), ytgt=tp->GetY( zK18Tgt );
      
      if( i<MaxHits ){
	event.nhBcOut[i]=nh;
	event.chisqrBcOut[i]=chisqr;
	event.u0BcOut[i]=u0;
	event.v0BcOut[i]=v0;
	event.x0BcOut[i]=x0;
	event.y0BcOut[i]=y0;
	event.xtgtBcOut[i]=xtgt;
	event.ytgtBcOut[i]=ytgt;
      }
      
      HF1( 1102, double(nh) );
      HF1( 1103, chisqr );
      HF1( 1104, x0 ); HF1( 1105, y0 );
      HF1( 1106, u0 ); HF1( 1107, v0 );
      HF2( 1108, x0, u0 ); HF2( 1109, y0, v0 );
      HF2( 1110, x0, y0 );
      HF1( 1111, xtgt ); HF1( 1112, ytgt );
      HF2( 1113, xtgt, ytgt );

      //evDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.DrawBcOutLocalTrack(tp);
      //       }
      //------------------------Cut
      if( MinTgtXBcOut<xtgt && xtgt<MaxTgtXBcOut && 
	  MinTgtYBcOut<ytgt && ytgt<MaxTgtYBcOut &&
	  chisqr<MaxChisqrBcOut ){
	++ntOk;
      }
    }
    if( ntOk<1 ) {
      //     if (FlagEvDisp) {
      //       std::cout << "rejected by BcOut" << std::endl;
      //       const EvDisp & evDisp = EvDisp::GetInstance();
      //       evDisp.get_command();
      //       evDisp.EndOfEvent();
      //     }
      return true;
    }
  }
  HF1( 1, 15. );

  //SdcOut Analysis
  DCAna->TrackSearchSdcOut();
  int ntSdcOut=DCAna->GetNtracksSdcOut();
  event.ntSdcOut=ntSdcOut;
  HF1( 1201, double(ntSdcOut) );
  //------------------------Cut
  if( ntSdcOut<1 || ntSdcOut>MaxNTrackSdcOut ) {
    //     if (FlagEvDisp) {
    //       std::cout << "rejected by SdcOut" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 16. );

  {
    int ntOk=0;
    for( int i=0; i<ntSdcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX( zTof ), y0=tp->GetY( zTof );
      
      if( i<MaxHits ){
	event.nhSdcOut[i]=nh;
	event.chisqrSdcOut[i]=chisqr;
	event.u0SdcOut[i]=u0;
	event.v0SdcOut[i]=v0;
	event.xtofSdcOut[i]=x0;
	event.ytofSdcOut[i]=y0;
      }
      
      HF1( 1202, double(nh) );
      HF1( 1203, chisqr );
      HF1( 1204, x0 ); HF1( 1205, y0 );
      HF1( 1206, u0 ); HF1( 1207, v0 );
      HF2( 1208, x0, u0 ); HF2( 1209, y0, v0 );
      HF2( 1210, x0, y0 );

      //evDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.DrawSdcOutLocalTrack(tp);
      //       }
      //------------------------Cut
      if( chisqr<MaxChisqrSdcOut ){
	++ntOk;
      }
    }
    if( ntOk<1 ) {
      //     if (FlagEvDisp) {
      //       std::cout << "rejected by SdcOut" << std::endl;
      //       const EvDisp & evDisp = EvDisp::GetInstance();
      //       evDisp.get_command();
      //       evDisp.EndOfEvent();
      //     }
      return true;
    }
  }
  HF1( 1, 17. );

  HF1( 1, 20. );

  //SKSTracking Analysis
  DCAna->TrackSearchSks();
  int ntSks=DCAna->GetNTracksSks();
  event.ntSks=ntSks;
  HF1( 3001, double(ntSks) );
  //------------------------Cut
  if( ntSks<1 || ntSks>MaxNTrackSks ) {
    //     if (FlagEvDisp) {
    //       std::cout << "rejected by SksTracking" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 21. );

  {
    int ntOk=0;
    for( int i=0; i<ntSks; ++i ){
      SksTrack *track=DCAna->GetSksTrack(i);
      if(!track) continue;
      int nh=track->GetNHits();
      double chisqr=track->chisqr();
      ThreeVector Ppos=track->PrimaryPosition();
      ThreeVector Pmom=track->PrimaryMomentum();
      double path=track->PathLengthToTOF();
      double xt=Ppos.x(), yt=Ppos.y();
      double p=Pmom.mag();
      double ut=Pmom.x()/p, vt=Pmom.y()/p;

      if(i<MaxHits){
	event.nhSks[i] = nh;
	event.chisqrSks[i] = chisqr;
	event.path[i] = path;
	event.pSks[i] = p;
	event.xts[i] = xt;
	event.yts[i] = yt;
	event.uts[i] = ut;
	event.vts[i] = vt;
      }

      HF1( 3002, double(nh) );
      HF1( 3003, chisqr );
      HF1( 3004, xt ); HF1( 3005, yt ); 
      HF1( 3006, ut ); HF1( 3007, vt );
      HF2( 3008, xt, ut ); HF2( 3009, yt, vt );
      HF2( 3010, xt, yt );
      HF1( 3011, p ); 
      HF1( 3012, path );

      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	double m2=MassSquare( p, path, clTof->CMeanTime()-time0+offset );
	// 	std::cout<<"Mom="<< p <<std::endl;
	// 	std::cout<<"Path="<< pathL <<std::endl;
	// 	std::cout<<"Time="<< clTof->CMeanTime()-time0 <<std::endl;
	std::cout<<"m2= "<< m2 <<std::endl;
	
	if( i<MaxHits ){
	  event.m2[i] = m2;
	}
	HF1( 3013, m2 );
      }
      //------------------------Cut
      if( chisqr<MaxChisqrSks ){
	++ntOk;
	HF1( 3102, double(nh) );
	HF1( 3103, chisqr );
	HF1( 3104, xt ); HF1( 3105, yt ); 
	HF1( 3106, ut ); HF1( 3107, vt );
	HF2( 3108, xt, ut ); HF2( 3109, yt, vt );
	HF2( 3110, xt, yt );
	HF1( 3111, p ); 
	HF1( 3112, path );

	for( int j=0; j<ncTof; ++j ){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if( !clTof || !clTof->GoodForAnalysis() ) continue;
	  double m2=MassSquare( p, path, clTof->CMeanTime()-time0+offset );
	  HF1( 3113, m2 );
	}
      }
      else {
	track->GoodForAnalysis( false );
      }
    }
    HF1( 3101, double(ntOk) );
    //------------------------Cut
    if( ntOk<1 ) {
      //       if (FlagEvDisp) {
      // 	std::cout << "rejected by Sks Chisq cut" << std::endl;
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.get_command();
      // 	evDisp.EndOfEvent();
      //       }
      return true;
    }
  }
  HF1( 1, 22. );

  HF1( 1, 30. );

  //BcIn Analysis
  DCAna->TrackSearchBcIn();
  int ntBcIn=DCAna->GetNtracksBcIn();
  event.ntBcIn=ntBcIn;
  HF1( 1301, double(ntBcIn) );
  //------------------------Cut
  if( ntBcIn<1 || ntBcIn>MaxNTrackBcIn ) {
    //     if (FlagEvDisp) {
    //       std::cout << "rejected by BcIn" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 31. );

  {
    int ntOk=0;
    for( int i=0; i<ntBcIn; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackBcIn(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(0.), y0=tp->GetY(0.);
      
      if(i<MaxHits){
	event.nhBcIn[i]=nh;
	event.chisqrBcIn[i]=chisqr;
	event.u0BcIn[i]=u0;
	event.v0BcIn[i]=v0;
	event.x0BcIn[i]=x0;
	event.y0BcIn[i]=y0;
      }
      
      HF1( 1302, double(nh) );
      HF1( 1303, chisqr );
      HF1( 1304, x0 ); HF1( 1305, y0 );
      HF1( 1306, u0 ); HF1( 1307, v0 );
      HF2( 1308, x0, u0 ); HF2( 1309, y0, v0 );
      HF2( 1310, x0, y0 );

      //evDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.DrawBcInLocalTrack(tp);
      //       } 
      //------------------------Cut
      if( chisqr<MaxChisqrBcIn ){
	++ntOk;
      }
   }  
    HF1( 3101, double(ntOk) );
    if( ntOk<1 ) {
      //     if (FlagEvDisp) {
      //       std::cout << "rejected by BcIn" << std::endl;
      //       const EvDisp & evDisp = EvDisp::GetInstance();
      //       evDisp.get_command();
      //       evDisp.EndOfEvent();
      //     }
      return true;
    }
  }
  HF1( 1, 32. );

  // K18 Tracking
  DCAna->TrackSearchK18();
  int ntK18=DCAna->GetNTracksK18();
  event.ntK18=ntK18;
  HF1( 2201, double(ntK18) );
  //------------------------Cut  
  if( ntK18<1 || ntK18>MaxNTrackK18) {
    //    if (FlagEvDisp) {
    //       std::cout << "rejected by K6Tracking" << std::endl;
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 33. );

  {
    int ntOk=0;
    for( int i=0; i<ntK18; ++i ){
      K18Track *track=DCAna->GetK18Track(i);
      if(!track) continue;
      int nh=track->GetNHitsTotal();
      double chisqr=track->chisquare();
      double pk18=track->P();
      
      if (i<MaxTrack) {
	event.nhK18[i]=nh;
	event.chisqrK18[i]=chisqr;
	event.pK18[i]=pk18;
      }

      HF1( 2202, double(nh) );
      HF1( 2203, chisqr );
      HF1( 2204, pk18 );

      ++ntOk;
      //------------------------Cut
      if( chisqr<MaxChisqrK18 ){
	++ntOk;
      }
    }
    HF1( 2251, double(ntOk) );
    if( ntOk<1 ) {
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.get_command();
      // 	evDisp.EndOfEvent();
      //       }
      return true;
    }
  }
  HF1( 1, 34. );
  
  HF1( 1, 40. );
  
  std::vector <ThreeVector> PionPCont, PionXCont;
  std::vector <ThreeVector> KaonPCont, KaonXCont;
  
  //Kaon cut
  for( int itSks=0; itSks<ntSks; ++itSks ){
    SksTrack *track=DCAna->GetSksTrack(itSks);
    //if( !track || !track->GoodForAnalysis() ) continue;
    //     DCLocalTrack *trIn =track->GetLocalTrackIn();
    //     DCLocalTrack *trOut=track->GetLocalTrackOut();
    HodoCluster *clTof=0;//, *clLc=0;
    //     double xtof=trOut->GetX(zTof), ytof=trOut->GetY(zTof);
    //     double xlc =trOut->GetX(zLc),  ylc =trOut->GetY(zLc);
    //     double mindifTof=100., mindifLc=200.;
    for( int j=0; j<ncTof; ++j ){
      HodoCluster *cl=hodoAna->GetClusterTOF(j);
      //if( !cl || !cl->GoodForAnalysis() ) continue;
      //       double x=(cl->MeanSeg()-7.5)*70.;
      //       if( fabs(x-xtof)<mindifTof ){
      clTof=cl; //mindifTof=fabs(x-xtof);
      //    }
    }
    //     for( int j=0; j<ncLc; ++j ){
    //       HodoCluster *cl=hodoAna->GetClusterLc(j);
    //       if( !cl || !cl->GoodForAnalysis() ) continue;
    //       double x=(cl->MeanSeg()-7.5)*100.;
    //       if( fabs(x-xlc)<mindifLc ){
    // 	clLc=cl; mindifLc=fabs(x-xlc);
    //       }
    //     }
    //    if( !clTof || !clLc ) continue;
    //if( !clTof ) continue;
    int nh=track->GetNHits();
    double chisqr=track->chisqr();
    ThreeVector Ppos=track->PrimaryPosition();
    ThreeVector Pmom=track->PrimaryMomentum();
    double pathL=track->PathLengthToTOF();
    double xt=Ppos.x(), yt=Ppos.y();
    double p=Pmom.mag();
    double ut=Pmom.x()/p, vt=Pmom.y()/p;
    double m2=MassSquare( p, pathL, clTof->CMeanTime()-time0+offset );
    HF1( 3202, double(nh) );
    HF1( 3203, chisqr );
    HF1( 3204, xt ); HF1( 3205, yt ); 
    HF1( 3206, ut ); HF1( 3207, vt );
    HF2( 3208, xt, ut ); HF2( 3209, yt, vt );
    HF2( 3210, xt, yt );
    HF1( 3211, p ); 
    HF1( 3212, pathL );
    HF1( 3213, m2 );

    //     double xTof=(clTof->MeanSeg()-7.5)*70.;
    //     double yTof=(clTof->TimeDif())*800./12.;
    //     double xLc=(clLc->MeanSeg()-7.5)*100.;
    //     double yLc=(clLc->TimeDif())*960./12.;
    //     HF2( 4011, clTof->MeanSeg()-0.5, xtof );
    //     HF2( 4012, clLc->MeanSeg()-0.5, xlc );
    //     HF2( 4013, clTof->TimeDif(), ytof );
    //     HF2( 4014, clLc->TimeDif(),  ylc );
    //     HF1( 4015, xtof-xTof ); HF1( 4017, ytof-yTof );
    //     HF1( 4016, xlc-xLc   ); HF1( 4018, ylc-yLc   );
    //     HF2( 4019, xtof-xTof, ytof-yTof );
    //     HF2( 4020, xlc-xLc,   ylc-yLc   );
    //------------------------Cut    
    if( m2>MinMassSquare && m2<MaxMassSquare ){
      //       double ttof=clTof->CMeanTime()-time0,
      // 	tlc=clLc->CMeanTime()-time0;
      //       HF1( 322, clTof->ClusterSize() );
      //       HF1( 323, clTof->MeanSeg()-0.5 );
      //       HF1( 324, ttof );
      //       HF1( 325, clTof->DeltaE() );
      //       HF1( 422, clLc->ClusterSize() );
      //       HF1( 423, clLc->MeanSeg()-0.5 );
      //       HF1( 424, tlc );
      //       HF1( 425, clLc->DeltaE() );
      //       HF2( 521, clTof->MeanSeg()-0.5, clLc->MeanSeg()-0.5 );
      //       HF2( 522, ttof, tlc );
      //       HF1( 523, tlc-ttof );
      
      //       double u0in=trIn->GetU0();
      //       HF2( 4001, u0in, ttof ); HF2( 4003, u0in, ttof+12.5*u0in );
      //       HF2( 4002, u0in, tlc  ); HF2( 4004, u0in, tlc+12.5*u0in );

      HF1( 4202, double(nh) );
      HF1( 4203, chisqr );
      HF1( 4204, xt ); HF1( 4205, yt ); 
      HF1( 4206, ut ); HF1( 4207, vt );
      HF2( 4208, xt, ut ); HF2( 4209, yt, vt );
      HF2( 4210, xt, yt );
      HF1( 4211, p ); 
      HF1( 4212, pathL );
      KaonPCont.push_back(Pmom); KaonXCont.push_back(Ppos);
    }
  }
  if( KaonPCont.empty() ) {
    //     if (FlagEvDisp) {
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    //return true;
  }

  HF1( 1, 41. );

  //Pion cut  
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    K18Track *track=DCAna->GetK18Track(itK18);
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trIn=track->TrackIn(), *trOut=track->TrackOut();
    if( !trIn || !trOut ) continue;
    double p=track->P(), x=track->Xtgt(), y=track->Ytgt();
    double u=track->Utgt(), v=track->Vtgt();
    double pt=p/sqrt(1.+u*u+v*v);
    ThreeVector Pos( x+x_off, y+y_off, 0. );
    ThreeVector Mom( pt*(u_off), pt*(v_off), pt );
    double xi=trIn->GetX0(), yi=trIn->GetY0();
    double ui=trIn->GetU0(), vi=trIn->GetV0();
    double xo=trOut->GetX0(), yo=trOut->GetY0();
    double chisqr=track->chisquare();
    //     double xbh2=(clBH2Time0->MeanSeg()-2.5)*15.+7.5;
    //     double xbh1dc=trIn->GetX(-700.), xbh2dc=trOut->GetX(680.);
    //     double mindif=2*MaxXDifBh1;
    //     HodoCluster *clBh1=0;
    //     for( int j=0; j<ncBh1; ++j ){
    //       HodoCluster *cl=hodoAna->GetClusterBH1(j);
    //       if( !cl || !cl->GoodForAnalysis() ) continue;
    //       double xb=(cl->MeanSeg()-4.)*50./3.;
    //       if( fabs(xb-xbh1dc)<mindif ) {
    // 	mindif=fabs(xb-xbh1dc); clBh1=cl;
    //       }
    //     }
    //     if(!clBh1) continue;
    //     double xbh1=(clBh1->MeanSeg()-4.)*50./3.;
   
    HF1( 4102, double(track->GetNHitsTotal()) );
    HF1( 4103, chisqr ); HF1( 4104, p ); 
    HF1( 4105, x ); HF1( 4106, y ); 
    HF1( 4107, xo ); HF1( 4108, yo ); HF1( 4109, u ); HF1( 4110, v );
    HF1( 4111, xi ); HF1( 4112, yi ); HF1( 4113, ui ); HF1( 4114, vi );

    //     HF1( 4121, clBH2Time0->MeanSeg()-0.5 );
    //     HF1( 4122, clBh1->MeanSeg()-0.5 );
    //     HF1( 4123, xbh2dc-xbh2 ); HF1( 4124, xbh1dc-xbh1 );
    //     HF1( 4125, clBH2Time0->DeltaE() ); HF1( 4126, clBh1->DeltaE() );
    //     HF1( 4127, time0-clBh1->CMeanTime() );

    PionPCont.push_back(Mom); PionXCont.push_back(Pos);
  }

  //MissingMass
  int nK=KaonPCont.size(), nPi=PionPCont.size();

  event.nK=nK;
  event.nPi=nPi;
  event.nPiK=nK*nPi;
  HF1( 4101, double(nPi));
  HF1( 4201, double(nK) );

  //------------------------Cut
  if( PionPCont.empty() ) {
    //     if (FlagEvDisp) {
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }
  HF1( 1, 42. );

  int npik=0;
  for( int ik=0; ik<nK; ++ik ){
    ThreeVector pk=KaonPCont[ik], xk=KaonXCont[ik];
    for( int ipi=0; ipi<nPi; ++ipi ){
      ThreeVector ppi=PionPCont[ipi], xpi=PionXCont[ipi];

      ThreeVector vert=VertexPoint( xpi, xk, ppi, pk );
      
      LorentzVector LvPi( ppi, sqrt(PionMass*PionMass+ppi.mag2()) );
      LorentzVector LvK( pk, sqrt(KaonMass*KaonMass+pk.mag2()) );
      LorentzVector LvC( 0., 0., 0., ProtonMass );
      LorentzVector LvCore( 0., 0., 0., 0. );
      LorentzVector LvRc=LvPi+LvC-LvK;
      double MisMass=LvRc.mag();//-LvC.mag();
      double BE=LvRc.mag()-( LvCore.mag()+PionMass );

     //double MisMass=LvRecoil.mag();
      double us=pk.x()/pk.z(), vs=pk.y()/pk.z();
      double ub=ppi.x()/ppi.z(), vb=ppi.y()/ppi.z();
      double cost=ppi*pk/(ppi.mag()*pk.mag());

      if (npik<MaxHits) {
	event.vtx[npik]=vert.x();
	event.vty[npik]=vert.y();
	event.vtz[npik]=vert.z();
	event.theta[npik]=acos(cost)*Rad2Deg;
	event.MissMass[npik]=MisMass;
	event.BE[npik]=BE;
	npik++;
      }
      //EvDisp
      //       if (FlagEvDisp) {
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	ThreeVector VertGlobal = 
      // 	  DCGeomMan::GetInstance().Local2GlobalPos(0,vert);
      // 	evDisp.DrawVertex(VertGlobal, cost);
      //       }
      
      HF1( 5001, vert.z() );
      //if( vert.z()<MinVertexZ || vert.z()>MaxVertexZ) continue;
      //------------------------Cut
      HF1( 1, 43. );

      HF1( 5002, MisMass );
      HF2( 5011, MisMass, us );
      HF2( 5012, MisMass, vs );
      HF2( 5013, MisMass, ub );
      HF2( 5014, MisMass, vb );
    }
  }  
  //   if (FlagEvDisp) {
  //     const EvDisp & evDisp = EvDisp::GetInstance();
  //     evDisp.get_command();
  //     evDisp.EndOfEvent();
  //   }

  HF1( 1, 50. );

  tree->Fill();

  return true;
}

void EventPiPiAna::InitializeEvent( void )
{

  //Trigger
  event.trigtype  = -1;

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }

  //Hodoscope
  for( int it=0; it<MaxHits; it++){
    event.bh1HitSeg[it] = -1;
    event.bh2HitSeg[it] = -1;
    event.tofHitSeg[it] = -1;
    event.lcHitSeg[it] = -1;

    event.tbh1Hit[it] = -1;
    event.tbh2Hit[it] = -1;
    event.ttofHit[it] = -1;
    event.tlcHit[it]  = -1;

    event.debh1Hit[it] = -1;
    event.debh2Hit[it] = -1;
    event.detofHit[it] = -1;
    event.delcHit[it]  = -1;

    event.btof[it]  = -1;
  }

  //DC
  event.ntBcIn   = -1;
  event.ntBcOut  = -1;
  event.ntSdcIn  = -1;
  event.ntSdcOut = -1;
  event.ntK18    = -1;
  event.ntSks    = -1;

  event.nPi      = -1;
  event.nK       = -1;
  event.nPiK     = -1;

  //Beam DC    
  for( int it=0; it<MaxHits; it++){
    event.nhBcIn[it]     = -1;
    event.chisqrBcIn[it] = -1.0;
    event.x0BcIn[it] = -9999.0;
    event.y0BcIn[it] = -9999.0;
    event.u0BcIn[it] = -9999.0;
    event.v0BcIn[it] = -9999.0;

    event.xbh1BcIn[it] = -9999.0;
    event.ybh1BcIn[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = -1;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.0;
    event.y0BcOut[it] = -9999.0;
    event.u0BcOut[it] = -9999.0;
    event.v0BcOut[it] = -9999.0;

    event.xtgtBcOut[it] = -9999.0;
    event.ytgtBcOut[it] = -9999.0;
    event.xbh2BcOut[it] = -9999.0;
    event.ybh2BcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = -1;
    event.chisqrK18[it] = -1.0;
    event.xtb[it] = -9999.0;
    event.ytb[it] = -9999.0;
    event.utb[it] = -9999.0;
    event.vtb[it] = -9999.0;

    event.pK18[it]   = -9999.0;
    event.thetab[it] = -9999.0;
  }

  //SKS DC
  for( int it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = -1;
    event.chisqrSdcIn[it] = -1.0;
    event.x0SdcIn[it] = -9999.0;
    event.y0SdcIn[it] = -9999.0;
    event.u0SdcIn[it] = -9999.0;
    event.v0SdcIn[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSdcOut[it]     = -1;
    event.chisqrSdcOut[it] = -1.0;
    event.u0SdcOut[it] = -9999.0;
    event.v0SdcOut[it] = -9999.0;

    event.xtofSdcOut[it] = -9999.0;
    event.ytofSdcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSks[it]     = -1;
    event.chisqrSks[it] = -1.0;
    event.xts[it] = -9999.0;
    event.yts[it] = -9999.0;
    event.uts[it] = -9999.0;
    event.vts[it] = -9999.0;

    event.pSks[it] = -9999.0;
    event.path[it] = -9999.0;
    event.m2[it]   = -9999.0;
    event.thetas[it] = -9999.0;
  }

  //Reaction
  for( int it=0; it<MaxHits; it++){
    event.vtx[it]      = -9999.0; 
    event.vty[it]      = -9999.0; 
    event.vtz[it]      = -9999.0; 
    event.theta[it]    = -9999.0; 
    event.MissMass[it] = -9999.0; 
    event.BE[it]       = -9999.0; 
  }
}

bool EventPiPiAna::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventPiPiAna;
}

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 113, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 114, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 115, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 201, "#Clusters BH1",  11, 0., 11. );
  HB1( 202, "ClusterSize BH1",11, 0., 11. );
  HB1( 203, "HitPat BH1", 11, 0., 11. );
  HB1( 204, "MeanTime BH1", 200, -10., 10. );
  HB1( 205, "Delta-E BH1", 200, -0.5, 4.5 );
  HB1( 206, "Beam ToF", 200, -10., 10. );

  HB1( 211, "#Clusters BH1 [pi]",  11, 0., 11. );
  HB1( 212, "ClusterSize BH1 [pi]",11, 0., 11. );
  HB1( 213, "HitPat BH1 [pi]", 11, 0., 11. );
  HB1( 214, "MeanTime BH1 [pi]", 200, -10., 10. );
  HB1( 215, "Delta-E BH1 [pi]", 200, -0.5, 4.5 );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -5., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -5., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 1001, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1002, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1003, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1004, "X0 SdcIn", 500, -200., 200. ); 
  HB1( 1005, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1006, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1007, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1102, "#Hits BcOut", 20, 0., 20. );
  HB1( 1103, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1104, "X0 BcOut", 500, -200., 200. ); 
  HB1( 1105, "Y0 BcOut", 500, -200., 200. );
  HB1( 1106, "U0 BcOut",  700, -0.35, 0.35 );
  HB1( 1107, "V0 BcOut",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200. );
  HB1( 1111, "Xtgt BcOut", 500, -200., 200. ); 
  HB1( 1112, "Ytgt BcOut", 500, -200., 200. );
  HB2( 1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

  HB1( 1201, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 1202, "#Hits SdcOut", 20, 0., 20. );
  HB1( 1203, "Chisqr SdcOut", 200, 0., 100. );
  HB1( 1204, "X0 SdcOut", 600, -1200., 1200. ); 
  HB1( 1205, "Y0 SdcOut", 600, -600., 600. );
  HB1( 1206, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 1207, "V0 SdcOut",  400, -0.20, 0.20 );
  HB2( 1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35 );
  HB2( 1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20 );
  HB2( 1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600. );

  HB1( 1301, "#Tracks BcIn", 10, 0., 10. );
  HB1( 1302, "#Hits BcIn", 20, 0., 20. );
  HB1( 1303, "Chisqr BcIn", 200, 0., 100. );
  HB1( 1304, "X0 BcIn", 500, -200., 200. ); 
  HB1( 1305, "Y0 BcIn", 500, -200., 200. );
  HB1( 1306, "U0 BcIn",  400, -0.35, 0.35 );
  HB1( 1307, "V0 BcIn",  200, -0.20, 0.20 );
  HB2( 1308, "X0%U0 BcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1309, "Y0%V0 BcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1310, "X0%Y0 BcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 2201, "#Tracks K18", 10, 0., 10. );
  HB1( 2202, "#Hits K18", 30, 0., 30. );
  HB1( 2203, "Chisqr K18", 500, 0., 50. );
  HB1( 2204, "P K18", 500, 0.2, 0.8 );
  HB1( 2251, "#Tracks K18 [Good]", 10, 0., 10. );

  HB1( 3001, "#Tracks Sks", 10, 0., 10. );
  HB1( 3002, "#Hits Sks", 30, 0., 30. );
  HB1( 3003, "Chisqr Sks", 500, 0., 500. );
  HB1( 3004, "Xtgt Sks", 500, -200., 200. );
  HB1( 3005, "Ytgt Sks", 500, -100., 100. );
  HB1( 3006, "Utgt Sks", 400, -0.35, 0.35 );
  HB1( 3007, "Vtgt Sks", 200, -0.20, 0.20 );
  HB2( 3008, "Xtgt%U Sks", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3009, "Ytgt%V Sks", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3010, "Xtgt%Ytgt Sks", 100, -200., 200., 100, -100., 100. );
  HB1( 3011, "P Sks", 200, 0.0, 1.0 );
  HB1( 3012, "PathLength Sks", 600, 3000., 6000. );
  HB1( 3013, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3101, "#Tracks Sks [Good]", 10, 0., 10. );
  HB1( 3102, "#Hits Sks [Good]", 30, 0., 30. );
  HB1( 3103, "Chisqr Sks [Good]", 500, 0., 500. );
  HB1( 3104, "Xtgt Sks [Good]", 500, -200., 200. );
  HB1( 3105, "Ytgt Sks [Good]", 500, -100., 100. );
  HB1( 3106, "Utgt Sks [Good]", 700, -0.35, 0.35 );
  HB1( 3107, "Vtgt Sks [Good]", 400, -0.20, 0.20 );
  HB2( 3108, "Xtgt%U Sks [Good]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3109, "Ytgt%V Sks [Good]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3110, "Xtgt%Ytgt Sks [Good]", 100, -200., 200., 100, -100., 100. );
  HB1( 3111, "P Sks [Good]", 200, 0.0, 1.0 );
  HB1( 3112, "PathLength Sks [Good]", 600, 3000., 6000. );
  HB1( 3113, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3201, "#Tracks Sks [Good2]", 10, 0., 10. );
  HB1( 3202, "#Hits Sks [Good2]", 30, 0., 30. );
  HB1( 3203, "Chisqr Sks [Good2]", 500, 0., 500. );
  HB1( 3204, "Xtgt Sks [Good2]", 500, -200., 200. );
  HB1( 3205, "Ytgt Sks [Good2]", 500, -100., 100. );
  HB1( 3206, "Utgt Sks [Good2]", 700, -0.35, 0.35 );
  HB1( 3207, "Vtgt Sks [Good2]", 400, -0.20, 0.20 );
  HB2( 3208, "Xtgt%U Sks [Good2]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3209, "Ytgt%V Sks [Good2]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3210, "Xtgt%Ytgt Sks [Good2]", 100, -200., 200., 100, -100., 100. );
  HB1( 3211, "P Sks [Good2]", 200, 0.0, 1.0 );
  HB1( 3212, "PathLength Sks [Good2]", 600, 3000., 6000. );
  HB1( 3213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 4101, "#Tracks K18 [SksP]", 10, 0., 10. );
  HB1( 4102, "#Hits K18 [SksP]", 30, 0., 30. );
  HB1( 4103, "Chisqr K18 [SksP]", 500, 0., 50. );
  HB1( 4104, "P K18 [SksP]", 200, 0.2, 0.8 );
  HB1( 4105, "Xtgt K18 [SksP]", 500, -200., 200. );
  HB1( 4106, "Ytgt K18 [SksP]", 500, -100., 100. );
  HB1( 4107, "Xout K18 [SksP]", 500, -200., 200. ); 
  HB1( 4108, "Yout K18 [SksP]", 500, -100., 100. );
  HB1( 4109, "Uout K18 [SksP]", 400, -0.35, 0.35 );
  HB1( 4110, "Vout K18 [SksP]", 200, -0.20, 0.20 );
  HB1( 4111, "Xin  K18 [SksP]", 500, -100., 100. ); 
  HB1( 4112, "Yin  K18 [SksP]", 500, -100., 100. );
  HB1( 4113, "Uin  K18 [SksP]", 400, -0.35, 0.35 );
  HB1( 4114, "Vin  K18 [SksP]", 200, -0.20, 0.20 );

  HB1( 4201, "#Tracks Sks [Proton]", 10, 0., 10. );
  HB1( 4202, "#Hits Sks [Proton]", 30, 0., 30. );
  HB1( 4203, "Chisqr Sks [Proton]", 500, 0., 500. );
  HB1( 4204, "Xtgt Sks [Proton]", 500, -200., 200. );
  HB1( 4205, "Ytgt Sks [Proton]", 500, -100., 100. );
  HB1( 4206, "Utgt Sks [Proton]", 700, -0.35, 0.35 );
  HB1( 4207, "Vtgt Sks [Proton]", 400, -0.20, 0.20 );
  HB2( 4208, "Xtgt%U Sks [Proton]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 4209, "Ytgt%V Sks [Proton]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 4210, "Xtgt%Ytgt Sks [Proton]", 100, -200., 200., 100, -100., 100. );
  HB1( 4211, "P Sks [Proton]", 200, 0.0, 1.0 );
  HB1( 4212, "PathLength Sks [Proton]", 600, 3000., 6000. );
  HB1( 4213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 5001, "Zvert [PiK]", 1000, -1000., 1000. ); 
  HB1( 5002, "MissingMass [PiK]", 500, 0.0, 0.5 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 0.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 0.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 0.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 0.50, 100, -0.10, 0.10 );
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");

  //Beam DC 
  tree->Branch("ntBcIn",    &event.ntBcIn,     "ntBcIn/I");
  tree->Branch("nhBcIn",     event.nhBcIn,     "nhBcIn[ntBcIn]/I");
  tree->Branch("chisqrBcIn", event.chisqrBcIn, "chisqrBcIn[ntBcIn]/D");
  tree->Branch("x0BcIn",     event.x0BcIn,     "x0BcIn[ntBcIn]/D");
  tree->Branch("y0BcIn",     event.y0BcIn,     "y0BcIn[ntBcIn]/D");
  tree->Branch("u0BcIn",     event.u0BcIn,     "u0BcIn[ntBcIn]/D");
  tree->Branch("v0BcIn",     event.v0BcIn,     "v0BcIn[ntBcIn]/D");
  tree->Branch("xbh1BcIn",   event.xbh1BcIn,   "xbh1BcIn[ntBcIn]/D");
  tree->Branch("ybh1BcIn",   event.ybh1BcIn,   "ybh1BcIn[ntBcIn]/D");

  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");
  tree->Branch("xtgtBcOut",  event.xtgtBcOut,   "xtgtBcOut[ntBcOut]/D");
  tree->Branch("ytgtBcOut",  event.ytgtBcOut,   "ytgtBcOut[ntBcOut]/D");
  tree->Branch("xbh2BcOut",  event.xbh2BcOut,   "xbh2BcOut[ntBcOut]/D");
  tree->Branch("ybh2BcOut",  event.ybh2BcOut,   "ybh2BcOut[ntBcOut]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");
  tree->Branch("xtb",         event.xtb,       "xtb[ntK18]/D");
  tree->Branch("ytb",         event.ytb,       "ytb[ntK18]/D");
  tree->Branch("utb",         event.utb,       "utb[ntK18]/D");
  tree->Branch("vtb",         event.vtb,       "vtb[ntK18]/D");
  tree->Branch("thetab",      event.thetab,    "thetab[ntK18]/D");

  //SKS
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrOut[ntSdcOut]/D");
  tree->Branch("xtofSdcOut",  event.xtofSdcOut,   "xtofSdcOut[ntSdcOut]/D");
  tree->Branch("ytofSdcOut",  event.ytofSdcOut,   "ytofSdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntSks",      &event.ntSks,     "ntSks/I");
  tree->Branch("nhSks",       event.nhSks,     "nhSks[ntSks]/I");
  tree->Branch("chisqrSks",   event.chisqrSks, "chisqrSks[ntSks]/D");
  tree->Branch("path",        event.path,      "path[ntSks]/D");
  tree->Branch("pSks",        event.pSks,      "pSks[ntSks]/D");
  tree->Branch("m2",          event.m2,        "m2[ntSks]/D");
  tree->Branch("xts",         event.xts,       "xts[ntSks]/D");
  tree->Branch("yts",         event.yts,       "yts[ntSks]/D");
  tree->Branch("uts",         event.uts,       "uts[ntSks]/D");
  tree->Branch("vts",         event.vts,       "vts[ntSks]/D");
  tree->Branch("thetas",      event.thetas,    "thetas[ntSks]/D");

  //Reaction
  tree->Branch("nPi",    &event.nPi,     "nPi/I");
  tree->Branch("nK",     &event.nK,      "nK/I");
  tree->Branch("nPiK",   &event.nPiK,    "nPiK/I");
  tree->Branch("vtx",     event.vtx,     "vtx[nPiK]/D");
  tree->Branch("vty",     event.vty,     "vty[nPiK]/D");
  tree->Branch("vtz",     event.vtz,     "vtz[nPiK]/D");
  tree->Branch("theta",   event.theta,    "theta[nPiK]/D");
  tree->Branch("MissMass",event.MissMass, "MissMass[nPiK]/D");
  tree->Branch("BE",      event.BE,       "BE[nPiK]/D");

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  if( HodoParamFileName_!="" )
    HodoParamManager_ = new HodoParamMan(HodoParamFileName_);
  if(HodoParamManager_) HodoParamManager_->Initialize();

  if( HodoPHCFileName_!="" )
    HodoPHCManager_ = new HodoPHCMan(HodoPHCFileName_);
  if( HodoPHCManager_ ) HodoPHCManager_->Initialize();

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();

  if( DCTdcCalibFileName_!="" )
    DCTdcCalibManager_ = new DCTdcCalibMan(DCTdcCalibFileName_);
  if(DCTdcCalibManager_) DCTdcCalibManager_->Initialize();

  if( DCDriftParamFileName_!="" )
    DCDriftParamManager_ = new DCDriftParamMan(DCDriftParamFileName_);
  if(DCDriftParamManager_) DCDriftParamManager_->Initialize();

  if( FieldMapFileName_!="" )
    FieldMan::GetInstance().Initialize( FieldMapFileName_ );

  if( K18MatrixFileName_!="" )
    K18Matrix_ = new K18TransMatrix(K18MatrixFileName_);
  if(K18Matrix_) K18Matrix_->Initialize();

  return true;
}
