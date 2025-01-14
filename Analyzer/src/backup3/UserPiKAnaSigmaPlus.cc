/*
  UserPiKAnaSigmaPlus.cc
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
#include "BH1Filter.hh"

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 2;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double offset = -3.0;

const double pB_offset   =  1.00114;
const double pS_offset   =  0.000;
const double pK18_offset =  0.00193403;
const double pSKS_offset = -0.00214165;
const double bh2_dE      =  0.0011;

const double x_off = -10.922;
const double y_off =  -0.6382;
const double u_off =   0.00183;
const double v_off =  -0.00219;

// const double x_off = 10.922;
// const double y_off =  0.6382;
// const double u_off = -0.00183;
// const double v_off =  0.00219;

VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventPiKAnaSigmaPlus 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventPiKAnaSigmaPlus();
  ~EventPiKAnaSigmaPlus();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventPiKAnaSigmaPlus::EventPiKAnaSigmaPlus()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventPiKAnaSigmaPlus::~EventPiKAnaSigmaPlus()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 50
#endif

//For Tree
struct Event{
  //Trigger
  int trigtype;
  int trigflag[NumOfMisc];

  //Hodoscope
  int nhitsBh2;
  int nhBH2;
  double BH2Seg[MaxHits];
  double tBH2[MaxHits];
  double deBH2[MaxHits];

  int nhBH1;
  double BH1Seg[MaxHits];
  double tBH1[MaxHits];
  double deBH1[MaxHits];
  double btof[MaxHits];

  int nhTof;
  double TofSeg[MaxHits];
  double tTof[MaxHits];
  double dtTof[MaxHits];
  double deTof[MaxHits];

  int nhLc;
  double LcSeg[MaxHits];
  double tLc[MaxHits];
  double dtLc[MaxHits];
  double deLc[MaxHits];

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
  double closeDist[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double BE[MaxHits];
  double theta_CM[MaxHits];
  double cost_CM[MaxHits];

  double xk[MaxHits];
  double yk[MaxHits];
  double uk[MaxHits];
  double vk[MaxHits];
  double xpi[MaxHits];
  double ypi[MaxHits];
  double upi[MaxHits];
  double vpi[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];
};
static Event event;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -1.0;
const double MaxBeamToF  =  1.0;

// TOF Cut Parameters
const double MinDeltaETof = 0.5; //for Kaon
const double MaxDeltaETof = 2.0; //for Kaon
const double MinTimeTof =  19.0;
const double MaxTimeTof =  30.0;

// LC Cut Parameters
const double MinDeltaELc = 0.00; //for Kaon
const double MaxDeltaELc = 3.50; //for Kaon
const double MinTimeLc =  20.0;
const double MaxTimeLc =  40.0;

// TOF-LC  
const double MinToFT2L =   0.0;
const double MaxToFT2L =  10.0;
const int MinToF2LSeg  =   0;
const int MaxToF2LSeg  =  10;

//SdcIn vs TOF&LC cut
//const double MinModTimeTof = 21.5;//For study
const double MinModTimeTof = 21.0;
const double MaxModTimeTof = 25.0;
const double MinModTimeLc  = 20.0;
const double MaxModTimeLc  = 30.0;

// TOF-SdcOut
const double MinXDifTof  =  -40.;
const double MaxXDifTof  =   70.;
const double MinYDifTof  =  -230.;
const double MaxYDifTof  =   200.;

// LC-SdcOut
const double MinXDifLc  = -80.;
const double MaxXDifLc  =  80.;
const double MinYDifLc  = -600.;
const double MaxYDifLc  =  400.;

//Multi cut of DCs
//For study
// const double MaxMultiHitBcIn   = 3.0;
// const double MaxMultiHitBcOut  = 10.;
// const double MaxMultiHitSdcIn  =  7.;
// const double MaxMultiHitSdcOut = 2.0;
const double MaxMultiHitBcIn   = 5.0;
const double MaxMultiHitBcOut  = 100.;
const double MaxMultiHitSdcIn  = 100.;
const double MaxMultiHitSdcOut = 3.0;

const double MaxChisqrBcIn   = 30.;
const double MaxChisqrBcOut  = 20.;
const double MaxChisqrSdcIn  = 20.;
const double MaxChisqrSdcOut = 20.;
const double MaxChisqrK18    = 100.;
const double MaxChisqrSks    = 100.;

const int MaxNTrackBcIn   = 30;
const int MaxNTrackBcOut  = 20;
const int MaxNTrackSdcIn  = 20;
const int MaxNTrackSdcOut = 20;
const int MaxNTrackK18    = 20;
const int MaxNTrackSks    = 20;

const double MinTgtXSdcIn  =  -100.;
const double MaxTgtXSdcIn  =   100.;
const double MinTgtYSdcIn  =  -50.;
const double MaxTgtYSdcIn  =   50.;

const double MinTgtXBcOut  =  -100.;
const double MaxTgtXBcOut  =   100.;
const double MinTgtYBcOut  =  -50.;
const double MaxTgtYBcOut  =   50.;

const double MinMassSquare =  0.10;
const double MaxMassSquare =  0.60;

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
const double SigmaPMass      = 1.18937;
const double SigmaMMass      = 1.19745;

bool EventPiKAnaSigmaPlus::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventPiKAnaSigmaPlus::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();
  
  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  HF1( 1, 0. );

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
      //if( seg==1 && T>0 ){ trig_type = seg; event.trigtype = seg;}
      //(pi, K)
      if( seg==2 && T>0 ){ trig_type = seg; event.trigtype = seg;}
      //(pi, pi)
      //if( seg==11 && T>0 ){ trig_type = seg; event.trigtype = seg;}

      event.trigflag[seg] = T;
    }
  }

  //------------------------Cut
  if (trig_type != PIK_TRIG )  return true;
  //if (trig_type != BEAM_TRIG )  return true;
  //if (trig_type != PIPI_TRIG )  return true;
  HF1( 1, 1. );

  //   std::cout << "PiK==" << trig_type <<std::endl;
  //   std::cout << "Flag==" << trig <<std::endl;

  ////////////////////////////Hodoscope
  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBH2=ncBh2;
  int nhBh2=hodoAna->GetNHitsBH2();
  event.nhitsBh2 = nhBh2;
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
    HF1( 103, clBH2Time0->MeanSeg()+1-0.5 );
    HF1( 104, mint );
    HF1( 105, clBH2Time0->DeltaE() );
    for( int i=1; i<ncBh2; ++i ){
      BH2Cluster *cl=hodoAna->GetClusterBH2(i);
      double t=cl->CMeanTime();
      double dEbh2 = cl->DeltaE();
      //------------------------Cut
      if( !(MinDeltaEBH2<dEbh2 && dEbh2<MaxDeltaEBH2) ) continue;

      HF1( 112, cl->ClusterSize() );
      HF1( 113, cl->MeanSeg()+1-0.5 );
      HF1( 114, t );
      HF1( 115, cl->DeltaE() );

      if(i<MaxHits){
	event.BH2Seg[i]=cl->MeanSeg()+1;
	event.tBH2[i]=t;
	event.deBH2[i]=cl->DeltaE();
      }
      if( fabs(t)<fabs(mint) ){
	clBH2Time0=cl;
	mint=t; time0=clBH2Time0->CTime0();
      }
    }
  }
  HF1( 122, clBH2Time0->ClusterSize() );
  HF1( 123, clBH2Time0->MeanSeg()+1-0.5 );
  HF1( 124, clBH2Time0->CMeanTime() );
  HF1( 125, clBH2Time0->DeltaE() );

  HF1( 1, 2. );
  
  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBH1=ncBh1;
  {
    int ncOk=0;
    HF1( 201, double(ncBh1) );
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      if( !(MinDeltaEBH1<dEbh1 && dEbh1<MaxDeltaEBH1) ) continue;

      HF1( 202, cl->ClusterSize() );
      HF1( 203, cl->MeanSeg()+1-0.5 );
      HF1( 204, cl->CMeanTime() );
      HF1( 205, cl->DeltaE() );
      HF1( 206, btof );

      if(i<MaxHits){
	event.BH1Seg[i]=cl->MeanSeg()+1;
	event.tBH1[i]=cl->CMeanTime();
	event.deBH1[i]=cl->DeltaE();
	event.btof[i] = btof;
      }
      //------------------------Cut
      if( btof>MinBeamToF && btof<MaxBeamToF ){
	++ncOk;
	HF1( 212, cl->ClusterSize() );
	HF1( 213, cl->MeanSeg()+1-0.5 );
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
  //       int seg=hit->SegmentId()+1;
  //       int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
  //       if (Tu>0 || Td>0) {
  // 	if (FlagEvDisp)	{
  // 	  const EvDisp & evDisp = EvDisp::GetInstance();	  
  // 	  evDisp.DrawHitHodoscope(TofId, seg+1, Tu, Td);
  // 	}
  //       }
  //     }
  //   }

  //////////////Tof Analysis
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  HF1( 301, double(ncTof) );
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0, dt=cl->TimeDif();
      HF1( 302, cl->ClusterSize() );
      HF1( 303, cl->MeanSeg()+1-0.5 );
      HF1( 304, t );
      HF1( 305, de );

      if(i<MaxHits){
	event.TofSeg[i]=cl->MeanSeg()+1;
	event.tTof[i]=t;
	event.dtTof[i]=dt;
	event.deTof[i]=de;
      }
      //------------------------Cut
      if( MinDeltaETof<de && de<MaxDeltaETof &&
	  MinTimeTof  <t  && t< MaxTimeTof && (cl->MeanSeg()+1)>12 ){
	++ncOk;
	HF1( 312, cl->ClusterSize() );
	HF1( 313, cl->MeanSeg()+1-0.5 );
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

  //////////////LC Analysis
  hodoAna->DecodeLCHits(rawData);
  int ncLc=hodoAna->GetNClustersLC();
  event.nhLc=ncLc;
  HF1( 401, double(ncLc) );
  {
    int ncOk=0;
    for( int i=0; i<ncLc; ++i ){
      HodoCluster *cl=hodoAna->GetClusterLC(i);
      double de=cl->DeltaE();
      double t=cl->CMeanTime()-time0, dt=cl->TimeDif();
      HF1( 402, cl->ClusterSize() );
      HF1( 403, cl->MeanSeg()+1-0.5 );
      HF1( 404, t );
      HF1( 405, de );

      if(i<MaxHits){
	event.LcSeg[i]=cl->MeanSeg()+1;
	event.tLc[i]=t;  
	event.dtLc[i]=dt;
	event.deLc[i]=de;
      }
      //------------------------Cut
      if( MinDeltaELc<de && de<MaxDeltaELc &&
	  MinTimeLc < t  && t <MaxTimeLc && (cl->MeanSeg()+1)>10 ){
	++ncOk;
	HF1( 412, cl->ClusterSize() );
	HF1( 413, cl->MeanSeg()+1-0.5 );
	HF1( 414, t );
	HF1( 415, de );
      }	
      else{
	cl->GoodForAnalysis( false );
      }
    }
    HF1( 411, double(ncOk) );
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
  HF1( 1, 5. );

  //////////////Tof-Lc
  {
    int nc=0;
    for( int iTof=0; iTof<ncTof; ++iTof ){
      HodoCluster *clTof=hodoAna->GetClusterTOF(iTof);
      if( !clTof || !clTof->GoodForAnalysis() ) continue;
      double ttof=clTof->CMeanTime()-time0,
	segTof=clTof->MeanSeg()+1;
      for( int iLc=0; iLc<ncLc; ++iLc ){
	HodoCluster *clLc=hodoAna->GetClusterLC(iLc);
	if( !clLc || !clLc->GoodForAnalysis() ) continue;
	double tlc=clLc->CMeanTime()-time0,
	  segLc=clLc->MeanSeg()+1;

	HF2( 501, segTof-0.5, segLc-0.5 );
	HF2( 502, ttof, tlc );
	HF1( 503, tlc-ttof );
	HF1( 504, segTof-segLc );
	//------------------------Cut
	if( MinToF2LSeg<=(segTof-segLc) && (segTof-segLc)<=MaxToF2LSeg &&
	    MinToFT2L<(tlc-ttof) && (tlc-ttof)<MaxToFT2L ){
	  ++nc;
	  HF2( 511, segTof-0.5, segLc-0.5 );
	  HF2( 512, ttof, tlc );
	  HF1( 513, tlc-ttof );
	  HF1( 514, segTof-segLc );
	}
	else{
	  clTof->GoodForAnalysis( false );
	  clLc->GoodForAnalysis( false );
	}
      }
    }
    HF1( 510, double(nc) );
    //------------------------Cut
    if( nc<1 ) {
      //       if (FlagEvDisp) {
      // 	std::cout << "rejected by Tof-Lc cut" << std::endl;
      // 	const EvDisp & evDisp = EvDisp::GetInstance();
      // 	evDisp.get_command();
      // 	evDisp.EndOfEvent();
      //       }
      return true;
    }
  }
  HF1( 1, 6. );

  //tree->Fill(); return true;
 
  HF1( 1, 10. );

  ////////////////////////////DC
  double multi_BcIn=0., multi_BcOut=0.;
  double multi_SdcIn=0., multi_SdcOut=0.;

  DCAna->DecodeRawHits( rawData ); 
  int IdTof = DCGeomMan::GetInstance().GetTofId(); 
  int IdLc  = DCGeomMan::GetInstance().GetLcId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 
  double zK18Tgt = DCGeomMan::GetInstance().GetLocalZ( IdK18Target );

  //////////////SDC3&4 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
    }
  }
  //------------------------Cut
  if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    return true;

  //////////////BC3&4 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  //------------------------Cut
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    return true;
  
  //////////////SDC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn  += double(nhIn);
    }
  }
  //------------------------Cut
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
    return true;
  
  HF1( 1, 11. );

  //////////////SdcIn Analysis
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
      
      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	double ttof=clTof->CMeanTime()-time0;
	//HF2( 1111, u0, ttof ); HF2( 1113, u0, ttof+12.5*u0 );
	//------------------------Cut
	if( MinModTimeTof< ttof+11.0*u0 && ttof+11.0*u0 <MaxModTimeTof )
	  condTof=true;
      }
      for( int j=0; j<ncLc; ++j ){
	HodoCluster *clLc=hodoAna->GetClusterLC(j);
	if( !clLc || !clLc->GoodForAnalysis() ) continue;
	double tlc=clLc->CMeanTime()-time0;
	//HF2( 1112, u0, tlc ); HF2( 1114, u0, tlc+12.5*u0 );
	//------------------------Cut
	//if( MinModTimeLc< tlc+10.0*u0 && tlc+10.0*u0 <MaxModTimeLc )
	condLc=true;
      }
      //------------------------Cut
      //if(  chisqr<MaxChisqrSdcIn ){
      if(  chisqr<MaxChisqrSdcIn && condTof && condLc ){
	++ntOk;
	// 	HF1( 1152, double(nh) );
	// 	HF1( 1153, chisqr );
	// 	HF1( 1154, x0 ); HF1( 1155, y0 );
	// 	HF1( 1156, u0 ); HF1( 1157, v0 );
	// 	HF2( 1158, x0, u0 ); HF2( 1159, y0, v0 );

	for( int j=0; j<ncTof; ++j ){
	  HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	  if( !clTof || !clTof->GoodForAnalysis() ) continue;
	  double ttof=clTof->CMeanTime()-time0;
	  //HF2( 1161, u0, ttof ); HF2( 1163, u0, ttof+12.5*u0 );
	}
	for( int j=0; j<ncLc; ++j ){
	  HodoCluster *clLc=hodoAna->GetClusterLC(j);
	  if( !clLc || !clLc->GoodForAnalysis() ) continue;
	  double tlc=clLc->CMeanTime()-time0;
	  //HF2( 1162, u0, tlc ); HF2( 1164, u0, tlc+12.5*u0 );
	}
	// For quick tracking
	//if(ntOk>0) tp->GoodForTracking( false );
      }
      else{
	tp->GoodForTracking( false );
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

  //////////////SdcOut Analysis
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
  HF1( 1, 14. );

  {
    int ntOk=0;
    for( int i=0; i<ntSdcOut; ++i ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(i);
      if(!tp) continue;
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX( zTof ), y0=tp->GetY( zTof );
      double xlc=tp->GetX( zLc ), ylc=tp->GetY( zLc );
    
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

      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
        if( !clTof || !clTof->GoodForAnalysis() ) continue;
        double x=(clTof->MeanSeg()-15.5)*70.;
        double y=(clTof->TimeDif())*50.;
	//HF1( TOFHid+41, x0-x ); HF1( TOFHid+42, y0-y );
	//------------------------Cut
	//if( MinXDifTof<x0-x && x0-x<MaxXDifTof )
	  // 	if( MinXDifTof<x0-x && x0-x<MaxXDifTof &&
	  //             MinYDifTof<y0-y && y0-y<MaxYDifTof )
          condTof=true;
      }
      for( int j=0; j<ncLc; ++j ){
        HodoCluster *clLc=hodoAna->GetClusterLC(j);
        if( !clLc || !clLc->GoodForAnalysis() ) continue;
        double x=(clLc->MeanSeg()-13.5)*100.;
        double y=(clLc->TimeDif())*50.;
	//HF1( LCHid+41, xlc-x ); HF1( LCHid+42, ylc-y );
	//------------------------Cut
	//if( MinXDifLc<xlc-x && xlc-x<MaxXDifLc )
// 	if( MinXDifLc<xlc-x && xlc-x<MaxXDifLc &&
//             MinYDifLc<ylc-y && ylc-y<MaxYDifLc )
          condLc=true;
      }
      if( condTof && condLc ){
	// 	HF1( 31, double(nh) ); HF1( 32, chisqr ); 
	// 	HF1( 34, x0 ); HF1( 35, y0 );
	// 	HF1( 36, u0 ); HF1( 37, v0 );
	// 	HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
	for( int ih=0; ih<nh; ++ih ){
	  DCLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer();
	  // 	  HF1( 33, layerId );
	}
        if( chisqr<MaxChisqrSdcOut ) ++ntOk;
	else 
	  tp->GoodForTracking( false );
      }
      else{
        tp->GoodForTracking( false );
      }
    }
    //    HF1( 30, double(ntOk) );
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
  HF1( 1, 15. );

  HF1( 1, 20. );

  //std::cout<<"***************************************"<<std::endl;
  //////////////SKSTracking Analysis
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

      double m2;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof || !clTof->GoodForAnalysis() ) continue;
	m2=MassSquare( p, path, clTof->CMeanTime()-time0+offset );
// 	std::cout<<"Mom="<< p <<std::endl;
// 	std::cout<<"Path="<< path <<std::endl;
// 	std::cout<<"Time="<< clTof->CMeanTime()-time0 <<std::endl;
// 	std::cout<<"m2= "<< m2 <<std::endl;
	
	if( i<MaxHits ){
	  event.m2[i] = m2;
	}
	HF1( 3013, m2 );
      }
      //------------------------Cut
      //if( chisqr<MaxChisqrSks ){
      if( MinMassSquare<m2 && m2<MaxMassSquare && chisqr<MaxChisqrSks ){
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
	  //std::cout<<"Time="<< clTof->CMeanTime()-time0+offset <<std::endl;
	  std::cout<<"******* m2= "<< m2 <<std::endl;
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

  //////////////BcOut Analysis
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
  HF1( 1, 31. );

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
      if( chisqr<MaxChisqrBcOut ){
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
  HF1( 1, 32. );

  //////////////BC1&2 number of hit in one layer
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcInHC(layer);
      int nhIn=contIn.size();
      multi_BcIn += double(nhIn);
    }
  }
  //------------------------Cut
  if( multi_BcIn/double(NumOfLayersBcIn) >MaxMultiHitBcIn )
    return true;

  HF1( 1, 33. );

  //////////////BcIn Analysis
  std::vector<std::vector<DCHitContainer> > bcInCandidates;
  BH1Filter& gFilter = BH1Filter::GetInstance();
  gFilter.Apply(*hodoAna, *DCAna, bcInCandidates);
  DCAna->TrackSearchBcIn(bcInCandidates);

  //DCAna->TrackSearchBcIn();
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
  HF1( 1, 34. );

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
  HF1( 1, 35. );

  //////////////K18 Tracking
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
  HF1( 1, 36. );

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
  HF1( 1, 37. );
  
  HF1( 1, 40. );

  std::vector <ThreeVector> PionPCont, PionXCont;
  std::vector <ThreeVector> KaonPCont, KaonXCont;
  
  //////////////Kaon cut
  for( int itSks=0; itSks<ntSks; ++itSks ){
    SksTrack *track=DCAna->GetSksTrack(itSks);
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trIn =track->GetLocalTrackIn();
    DCLocalTrack *trOut=track->GetLocalTrackOut();
    HodoCluster *clTof=0, *clLc=0;
    double xtof=trOut->GetX(zTof), ytof=trOut->GetY(zTof);
    double xlc =trOut->GetX(zLc),  ylc =trOut->GetY(zLc);
    double mindifTof=MinXDifTof, mindifLc=MinXDifLc;
    for( int j=0; j<ncTof; ++j ){
      HodoCluster *cl=hodoAna->GetClusterTOF(j);
      if( !cl || !cl->GoodForAnalysis() ) continue;
      double x=(cl->MeanSeg()-15.5)*70.;
      //if( fabs(x-xtof)<mindifTof ){
      clTof=cl; //mindifTof=fabs(x-xtof);
      //}
    }
    for( int j=0; j<ncLc; ++j ){
      HodoCluster *cl=hodoAna->GetClusterLC(j);
      if( !cl || !cl->GoodForAnalysis() ) continue;
      double x=(cl->MeanSeg()-13.5)*100.;
      //if( fabs(x-xlc)<mindifLc ){
	clLc=cl; //mindifLc=fabs(x-xlc);
	//}
    }
    if( !clTof || !clLc ) continue;
    int nh=track->GetNHits();
    double chisqr=track->chisqr();
    ThreeVector Ppos=track->PrimaryPosition();
    ThreeVector Pmom=track->PrimaryMomentum();

    //Calibration
    double p0 = Pmom.mag()+pSKS_offset+pS_offset;
    //double p0 = Pmom.mag()/(1.-(pSKS_offset/Pmom.mag()));

    double u0 = Pmom.x()/p0, v0 = Pmom.y()/p0;
    double pt = p0/sqrt(1.+u0*u0+v0*v0);

    ThreeVector PposCorr( Ppos.x()+x_off, Ppos.y()+y_off, Ppos.z() );
    ThreeVector PmomCorr( pt*(u0+u_off), pt*(v0+v_off), pt);

    double pathL=track->PathLengthToTOF();
    double xt=PposCorr.x(), yt=PposCorr.y();
    double p=PmomCorr.mag();
    double ut=PmomCorr.x()/p, vt=PmomCorr.y()/p;
    double m2=MassSquare( p0, pathL, clTof->CMeanTime()-time0+offset );
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
    if( MinMassSquare<m2 && m2<MaxMassSquare ){
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
      KaonPCont.push_back(PmomCorr); KaonXCont.push_back(PposCorr);
    }
  }
  if( KaonPCont.empty() ) {
    //     if (FlagEvDisp) {
    //       const EvDisp & evDisp = EvDisp::GetInstance();
    //       evDisp.get_command();
    //       evDisp.EndOfEvent();
    //     }
    return true;
  }

  HF1( 1, 41. );

  //Pion cut  
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    K18Track *track=DCAna->GetK18Track(itK18);
    if( !track || !track->GoodForAnalysis() ) continue;
    DCLocalTrack *trIn=track->TrackIn(), *trOut=track->TrackOut();
    if( !trIn || !trOut ) continue;

    //Calibration
    double p_offset;
    if( nhBh2==2 ) p_offset = pK18_offset - bh2_dE;
    else p_offset = pK18_offset;

    double p=track->P()*pB_offset+p_offset, x=track->Xtgt(), y=track->Ytgt();
    double u=track->Utgt(), v=track->Vtgt();
    //     double p=track->P(), x=track->Xtgt()+x_off, y=track->Ytgt()+y_off;
    //     double u=track->Utgt()+u_off, v=track->Vtgt()+v_off;

    double pt=p/sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    double xi=trIn->GetX0(), yi=trIn->GetY0();
    double ui=trIn->GetU0(), vi=trIn->GetV0();
    double xo=trOut->GetX0(), yo=trOut->GetY0();
    double chisqr=track->chisquare();
    //     double xbh2=(clBH2Time0->MeanSeg()-2.5)*15.+7.5;
    //     double xbh1dc=trIn->GetX(-700.), xbh2dc=trOut->GetX(680.);
    //     double mindif=2*MaxXDifBh1;
    HodoCluster *clBh1=0;
    for( int j=0; j<ncBh1; ++j ){
      HodoCluster *cl=hodoAna->GetClusterBH1(j);
      //if( !cl || !cl->GoodForAnalysis() ) continue;
    //       double xb=(cl->MeanSeg()-4.)*50./3.;
    //       if( fabs(xb-xbh1dc)<mindif ) {
    // 	mindif=fabs(xb-xbh1dc);
      clBh1=cl;
    }
    //     }
    if(!clBh1) continue;
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

  //MissingMass
  int nK=KaonPCont.size(), nPi=PionPCont.size();

  event.nK=nK;
  event.nPi=nPi;
  event.nPiK=nK*nPi;
  HF1( 4101, double(nPi));
  HF1( 4201, double(nK) );

  int npik=0;
  for( int ik=0; ik<nK; ++ik ){
    ThreeVector pk=KaonPCont[ik], xk=KaonXCont[ik];
    for( int ipi=0; ipi<nPi; ++ipi ){
      ThreeVector ppi=PionPCont[ipi], xpi=PionXCont[ipi];

      ThreeVector vert=VertexPoint( xpi, xk, ppi, pk );
      double closedist=closeDist( xpi, xk, ppi, pk );

      double us=pk.x()/pk.z(), vs=pk.y()/pk.z();
      double ub=ppi.x()/ppi.z(), vb=ppi.y()/ppi.z();
      double cost=ppi*pk/(ppi.mag()*pk.mag());

      //For correction @ QM02b for narrow region, 4th pol for V
      double pu0 =  0.0;
      double pu1 =  0.00268264;
      double pu2 =  0.00958551;
      double pu3 = -0.0627957;
      double pu4 =  0.40378;

      double pv0 =  0.000034491;
      double pv1 = -0.000914295;
      double pv2 = -0.000035015;
      double pv3 =  0.117698;
      double pv4 =  4.57913;

      //       //For correction @ QM02 for narrow region, 4th pol for V
      //       double pu0 =  0.0;
      //       double pu1 =  0.00245557;
      //       double pu2 =  0.00897695;
      //       double pu3 = -0.113089;
      //       double pu4 =  0.432664;
      
      //       double pv0 =  0.000077;
      //       double pv1 = -0.00238278;
      //       double pv2 =  0.0537012;
      //       double pv3 =  0.379753;
      //       double pv4 =-11.0171;

      double pk0= pk.mag();
      double pCorr = pk0 + (( pu1*us + pu2*us*us + pu3*us*us*us + pu4*us*us*us*us )
			    + (  pv1*vs + pv2*vs*vs + pv3*vs*vs*vs + pv4*vs*vs*vs*vs )) + (pu0+pv0);

      ThreeVector pkCorr(pCorr*pk.x()/pk.mag(),
			 pCorr*pk.y()/pk.mag(),
			 pCorr*pk.z()/pk.mag());

      //ThreeVector vert_eLoss = ( vert.x(), vert.y(), vert.z()+10.0 );

      ThreeVector ppiCorrDE = CorrElossIn( ppi, xpi, vert, PionMass);
      ThreeVector pkCorrDE  = CorrElossOut( pkCorr, xk, vert, KaonMass);

      LorentzVector LvPi( ppi, sqrt(PionMass*PionMass+ppi.mag2()) );
      LorentzVector LvPiCorrDE( ppiCorrDE, sqrt(PionMass*PionMass+ppiCorrDE.mag2()) );

      LorentzVector LvK( pk, sqrt(KaonMass*KaonMass+pk.mag2()) );
      LorentzVector LvKCorr( pkCorr, sqrt(KaonMass*KaonMass+pkCorr.mag2()) );
      LorentzVector LvKCorrDE( pkCorrDE, sqrt(KaonMass*KaonMass+pkCorrDE.mag2()) );

      LorentzVector LvC( 0., 0., 0., ProtonMass );
      LorentzVector LvCore( 0., 0., 0., 0. );

      LorentzVector LvRc = LvPi+LvC-LvK;
      LorentzVector LvRcCorr = LvPi+LvC-LvKCorr;
      LorentzVector LvRcCorrDE = LvPiCorrDE+LvC-LvKCorrDE;
      double MisMass = LvRc.mag();//-LvC.mag();
      double MisMassCorr = LvRcCorr.mag();//-LvC.mag();
      double MisMassCorrDE = LvRcCorrDE.mag();//-LvC.mag();
      double BE=LvRc.mag()-( LvCore.mag()+KaonMass );
      double BECorr=LvRcCorr.mag()-( LvCore.mag()+KaonMass );
      double BECorrDE=LvRcCorrDE.mag()-( LvCore.mag()+KaonMass );

      std::cout<<"******* Missing Mass= "<< MisMass 
	       << " (" << MisMassCorr << ")" 
	       << " (" << MisMassCorrDE << ")" <<std::endl;

      //Primary frame
      LorentzVector PrimaryLv =  LvPi+LvC;
      double TotalEnergyCM = PrimaryLv.mag();
      ThreeVector beta( PrimaryLv.vect()/PrimaryLv.e() );
      
      //CM
      double TotalMomCM 
	= 0.5*sqrt(( TotalEnergyCM*TotalEnergyCM
		     -( KaonMass+SigmaPMass )*( KaonMass+SigmaPMass ))
		   *( TotalEnergyCM*TotalEnergyCM
		      -( KaonMass-SigmaPMass )*( KaonMass-SigmaPMass )))/TotalEnergyCM;
      
      double costLab = cost;
      double cottLab = costLab/sqrt(1.-costLab*costLab);
      double bt=beta.mag(), gamma=1./sqrt(1.-bt*bt);
      double gbep=gamma*bt*sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
      double a  = gamma*gamma+cottLab*cottLab;
      double bp = gamma*gbep;
      double c  = gbep*gbep-cottLab*cottLab;
      double dd=bp*bp-a*c;
      
      if( dd<0. ){
	std::cerr << "dd<0." << std::endl;
	dd=0.;
	//exit(-1);
      }
      
      double costCM=(sqrt(dd)-bp)/a;
      if( costCM>1. || costCM<-1. ){
	std::cerr << "costCM>1. || costCM<-1." << std::endl;
	//exit(-1);
	costCM=-1.;
      }
      double sintCM=sqrt(1.-costCM*costCM);
      double KaonMom = TotalMomCM*sintCM/sqrt(1.-costLab*costLab);

      if (npik<MaxHits) {
	event.vtx[npik]=vert.x();
	event.vty[npik]=vert.y();
	event.vtz[npik]=vert.z();
	event.closeDist[npik]=closedist;
	event.theta[npik]=acos(cost)*Rad2Deg;
	event.theta_CM[npik] =acos(costCM)*Rad2Deg;
	event.cost_CM[npik] =costCM;

	event.MissMass[npik]=MisMass;
	event.MissMassCorr[npik] = MisMassCorr;
	event.MissMassCorrDE[npik] = MisMassCorrDE;

	event.BE[npik]=BE;

	event.xk[npik] = xk.x();
	event.yk[npik] = xk.y();
	event.uk[npik] = us;
	event.vk[npik] = vs;

	event.xpi[npik] = xpi.x();
	event.ypi[npik] = xpi.y();
	event.upi[npik] = ub;
	event.vpi[npik] = vb;
	event.pOrg[npik] = pk0;
	event.pCalc[npik] = KaonMom;
	event.pCorr[npik] = pCorr;
	event.pCorrDE[npik] = pkCorrDE.mag();
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

  //if( npik==0 ) return true;

  HF1( 1, 50. );

  tree->Fill();

  //Final Hodoscope histograms
  for( int i=1; i<ncBh2; ++i ){
    BH2Cluster *cl=hodoAna->GetClusterBH2(i);
    double t=cl->CMeanTime();
    double dEbh2 = cl->DeltaE();
    if( !(MinDeltaEBH2<dEbh2 && dEbh2<MaxDeltaEBH2) ) continue;

    HF1( 152, cl->ClusterSize() );
    HF1( 153, cl->MeanSeg()+1-0.5 );
    HF1( 154, t );
    HF1( 155, cl->DeltaE() );
  }

  for( int i=0; i<ncBh1; ++i ){
    HodoCluster *cl=hodoAna->GetClusterBH1(i);
    double btof= (cl->CMeanTime())-time0;
    if( !cl || !cl->GoodForAnalysis() ) continue;

    HF1( 252, cl->ClusterSize() );
    HF1( 253, cl->MeanSeg()+1-0.5 );
    HF1( 254, cl->CMeanTime() );
    HF1( 255, cl->DeltaE() );
    HF1( 256, btof );
  }

  for( int i=0; i<ncTof; ++i ){
    HodoCluster *cl=hodoAna->GetClusterTOF(i);
    double de=cl->DeltaE();
    double t=cl->CMeanTime()-time0;
    if( !cl || !cl->GoodForAnalysis() ) continue;

    HF1( 352, cl->ClusterSize() );
    HF1( 353, cl->MeanSeg()+1-0.5 );
    HF1( 354, t );
    HF1( 355, de );
  }
  
  for( int i=0; i<ncLc; ++i ){
    HodoCluster *cl=hodoAna->GetClusterLC(i);
    double de=cl->DeltaE();
    double t=cl->CMeanTime()-time0;
    if( !cl || !cl->GoodForAnalysis() ) continue;

    HF1( 452, cl->ClusterSize() );
    HF1( 453, cl->MeanSeg()+1-0.5 );
    HF1( 454, t );
    HF1( 455, de );
  }

  for( int iTof=0; iTof<ncTof; ++iTof ){
    HodoCluster *clTof=hodoAna->GetClusterTOF(iTof);
    if( !clTof || !clTof->GoodForAnalysis() ) continue;
    double ttof=clTof->CMeanTime()-time0,
      segTof=clTof->MeanSeg()+1;
    for( int iLc=0; iLc<ncLc; ++iLc ){
      HodoCluster *clLc=hodoAna->GetClusterLC(iLc);
      if( !clLc || !clLc->GoodForAnalysis() ) continue;
      double tlc=clLc->CMeanTime()-time0,
	segLc=clLc->MeanSeg()+1;
      
      HF2( 551, segTof-0.5, segLc-0.5 );
      HF2( 552, ttof, tlc );
      HF1( 553, tlc-ttof );
      HF1( 554, segTof-segLc );
    }
  }

  return true;
}

void EventPiKAnaSigmaPlus::InitializeEvent( void )
{

  //Trigger
  event.trigtype  = -1;

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }

  //Hodoscope
  event.nhitsBh2  = -1;
  event.nhBH2  = -1;
  event.nhBH1  = -1;
  event.nhTof  = -1;
  event.nhLc   = -1;

  for( int it=0; it<MaxHits; it++){
    event.BH2Seg[it] = -1;
    event.tBH2[it] = -9999.0;
    event.deBH2[it] = -9999.0;

    event.BH1Seg[it] = -1;
    event.tBH1[it] = -9999.0;
    event.deBH1[it] = -9999.0;
    event.btof[it]  = -1;

    event.TofSeg[it] = -1;
    event.tTof[it] = -9999.0;
    event.dtTof[it] = -9999.0;
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;   
    event.dtLc[it] = -9999.0;
    event.deLc[it] = -9999.0;
  }

  //DC
  event.ntBcIn   = -1;
  event.ntBcOut  = -1;
  event.ntSdcIn  = -1;
  event.ntSdcOut = -1;
  event.ntK18    = -1;
  event.ntSks    = -1;

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
  event.nPi      = -1;
  event.nK       = -1;
  event.nPiK     = -1;

  for( int it=0; it<MaxHits; it++){
    event.vtx[it]       = -9999.0; 
    event.vty[it]       = -9999.0; 
    event.vtz[it]       = -9999.0; 
    event.closeDist[it] = -9999.0; 
    event.theta[it]     = -9999.0; 
    event.theta_CM[it]  = -9999.0; 
    event.cost_CM[it]   = -9999.0; 
    event.MissMass[it]  = -9999.0;
    event.MissMassCorr[it]  = -9999.0; 
    event.MissMassCorrDE[it]  = -9999.0;
    event.BE[it]        = -9999.0; 

    event.xk[it] = -9999.0;
    event.yk[it] = -9999.0;
    event.uk[it] = -9999.0;
    event.vk[it] = -9999.0;
    event.xpi[it] = -9999.0;
    event.ypi[it] = -9999.0;
    event.upi[it] = -9999.0;
    event.vpi[it] = -9999.0;

    event.pOrg[it] = -9999.0;
    event.pCalc[it] = -9999.0;
    event.pCorr[it] = -9999.0;
    event.pCorrDE[it] = -9999.0;
  }
}

bool EventPiKAnaSigmaPlus::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventPiKAnaSigmaPlus;
}

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2", 7, 0., 7. );
  HB1( 113, "HitPat BH2", 10, 0., 10. );
  HB1( 114, "MeanTime BH2", 200, -10., 10. );
  HB1( 115, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 122, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 123, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 124, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 125, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 152, "ClusterSize BH2 [piK]", 7, 0., 7. );
  HB1( 153, "HitPat BH2 [piK]", 10, 0., 10. );
  HB1( 154, "MeanTime BH2 [piK]", 200, -10., 10. );
  HB1( 155, "Delta-E BH2 [piK]", 200, -0.5, 4.5 );

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

  HB1( 252, "ClusterSize BH1 [piK]",11, 0., 11. );
  HB1( 253, "HitPat BH1 [piK]", 11, 0., 11. );
  HB1( 254, "MeanTime BH1 [piK]", 200, -10., 10. );
  HB1( 255, "Delta-E BH1 [piK]", 200, -0.5, 4.5 );
  HB1( 256, "Beam ToF [piK]", 200, -10., 10. );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -50., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -50., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 352, "ClusterSize Tof [piK]",32, 0., 32. );
  HB1( 353, "HitPat Tof [piK]", 32, 0., 32. );
  HB1( 354, "TimeOfFlight Tof [piK]", 500, -50., 100. );
  HB1( 355, "Delta-E Tof [piK]", 200, -0.5, 4.5 );

  HB1( 401, "#Clusters Lc",  28, 0., 28. );
  HB1( 402, "ClusterSize Lc",28, 0., 28. );
  HB1( 403, "HitPat Lc", 28, 0., 28. );
  HB1( 404, "TimeOfFlight Lc", 500, -50., 200. );
  HB1( 405, "Delta-E Lc", 200, -0.5, 4.5 );

  HB1( 411, "#Clusters Lc [Good]",  28, 0., 28. );
  HB1( 412, "ClusterSize Lc [Good]",28, 0., 28. );
  HB1( 413, "HitPat Lc [Good]", 28, 0., 28. );
  HB1( 414, "TimeOfFlight Lc [Good]", 500, -50., 200. );
  HB1( 415, "Delta-E Lc [Good]", 200, -0.5, 4.5 );

  HB1( 452, "ClusterSize Lc [piK]",28, 0., 28. );
  HB1( 453, "HitPat Lc [piK]", 28, 0., 28. );
  HB1( 454, "TimeOfFlight Lc [piK]", 500, -50., 200. );
  HB1( 455, "Delta-E Lc [piK]", 200, -0.5, 4.5 );

  HB2( 501, "SegLC%SegTOF", 32, 0., 32., 28, 0., 28. );
  HB2( 502, "TLc%TTof", 100,-50., 100., 100, -50., 200. );
  HB1( 503, "TLc-TTof", 500, -30., 100. );
  HB1( 504, "TofSeg-LCSeg", 25, -10., 15. );

  HB1( 510, "#GoodCombi Tof-LC", 20, 0., 20. ); 
  HB2( 511, "SegLC%SegTOF [Good]", 32, 0., 32., 28, 0., 28. );
  HB2( 512, "TLc%TTof [Good]", 100,-20., 100., 100, -20., 100. );
  HB1( 513, "TLc-TTof [Good]", 500, -5., 100. );
  HB1( 514, "TofSeg-LCSeg [Good]", 25, -10., 15. );

  HB2( 551, "SegLC%SegTOF [piK]", 32, 0., 32., 28, 0., 28. );
  HB2( 552, "TLc%TTof [piK]", 100,-20., 100., 100, -20., 100. );
  HB1( 553, "TLc-TTof [piK]", 500, -5., 100. );
  HB1( 554, "TofSeg-LCSeg [piK]", 25, -10., 15. );

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
  HB1( 5002, "MissingMass [PiK]", 1000, 0.0, 2.0 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 2.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 2.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 2.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 2.50, 100, -0.10, 0.10 );
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Trigger  
  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");

  //Hodoscope
  tree->Branch("nhitsBh2",&event.nhitsBh2, "nhitsBh2/I");
  tree->Branch("nhBH2",   &event.nhBH2,   "nhBH2/I");
  tree->Branch("BH2Seg",   event.BH2Seg,  "BH2Seg[nhBH2]/D");
  tree->Branch("tBH2",     event.tBH2,    "tBH2[nhBH2]/D");
  tree->Branch("deBH2",    event.deBH2,   "deBH2[nhBH2]/D");

  tree->Branch("nhBH1",   &event.nhBH1,   "nhBH1/I");
  tree->Branch("BH1Seg",   event.BH1Seg,  "BH1Seg[nhBH1]/D");
  tree->Branch("tBH1",     event.tBH1,    "tBH1[nhBH1]/D");
  tree->Branch("deBH1",    event.deBH1,   "deBH1[nhBH1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBH1]/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhLc",   &event.nhLc,   "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/D");
  tree->Branch("tLc",     event.tLc,    "tLc[nhLc]/D");
  tree->Branch("dtLc",    event.dtLc,   "dtLc[nhLc]/D");
  tree->Branch("deLc",    event.deLc,   "deLc[nhLc]/D");

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
  tree->Branch("closeDist",event.closeDist,"closeDist[nPiK]/D");
  tree->Branch("theta",   event.theta,    "theta[nPiK]/D");
  tree->Branch("MissMass",event.MissMass, "MissMass[nPiK]/D");
  tree->Branch("MissMassCorr",event.MissMassCorr, "MissMassCorr[nPiK]/D");
  tree->Branch("MissMassCorrDE",event.MissMassCorrDE, "MissMassCorrDE[nPiK]/D");
  tree->Branch("BE",      event.BE,       "BE[nPiK]/D");
  tree->Branch("theta_CM", event.theta_CM, "theta_CM[nPiK]/D");
  tree->Branch("cost_CM",  event.cost_CM,  "cost_CM[nPiK]/D");

  tree->Branch("xk",         event.xk,       "xk[nPiK]/D");
  tree->Branch("yk",         event.yk,       "yk[nPiK]/D");
  tree->Branch("uk",         event.uk,       "uk[nPiK]/D");
  tree->Branch("vk",         event.vk,       "vk[nPiK]/D");
  tree->Branch("xpi",        event.xpi,      "xpi[nPiK]/D");
  tree->Branch("ypi",        event.ypi,      "ypi[nPiK]/D");
  tree->Branch("upi",        event.upi,      "upi[nPiK]/D");
  tree->Branch("vpi",        event.vpi,      "vpi[nPiK]/D");

  tree->Branch("pOrg",       event.pOrg,      "pOrg[nPiK]/D");
  tree->Branch("pCalc",      event.pCalc,     "pCalc[nPiK]/D");
  tree->Branch("pCorr",      event.pCorr,     "pCorr[nPiK]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,   "pCorrDE[nPiK]/D");

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

  if ( bh1FilterFileName_!="" )
    BH1Filter::GetInstance().Initialize(bh1FilterFileName_);

  return true;
}
