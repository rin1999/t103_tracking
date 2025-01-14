/*
  UserSksHodoscopeSigmaCut.cc
*/

#include "VEvent.hh"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "SksLib.hh"

const int BEAM_TRIG  = 1;
const int PIK_TRIG   = 2;
const int BTOF_TRIG  = 3;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;
const int BH2_TRIG   = 4;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double offset = -3.0;

static const int BH2Hid = 20000;
static const int TOFHid = 50000;
static const int ACHid  = 60000;
static const int LCHid  = 70000;

const double LightVel = 299.792458;  /* mm/nsec */
const double PionMass   = 0.139570;
const double KaonMass   = 0.493677;
const double ProtonMass = 0.9382720;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -1.0;
const double MaxBeamToF  =  1.0;

// DC multiplisity cut
const double MaxMultiHitSdcIn  = 100.0;
const double MaxMultiHitSdcOut =   3.0;

///////////////////////////////Tight
const double MinDeltaETof = 0.50; //for Kaon
const double MaxDeltaETof = 2.00; //for Kaon
const double MinTimeTof =  19.0;
const double MaxTimeTof =  30.0;

const double MinDeltaELc = 0.00; //for Kaon
const double MaxDeltaELc = 3.50; //for Kaon
const double MinTimeLc =  20.0;
const double MaxTimeLc =  40.0;

const double MinToFT2L =   0.0;
const double MaxToFT2L =  10.0;
const int MinToF2LSeg  =    0;
const int MaxToF2LSeg  =   10;

const double MinModTimeTof = 21.0;
const double MaxModTimeTof = 25.0;
const double MinModTimeLc  = 20.0;
const double MaxModTimeLc  = 30.0;

//Others
// TOF-SdcOut
const double MinXDifTof  =   -40.;
const double MaxXDifTof  =    70.;
const double MinYDifTof  =  -230.;
const double MaxYDifTof  =   200.;

// LC-SdcOut
const double MinXDifLc  =  -80.;
const double MaxXDifLc  =   80.;
const double MinYDifLc  = -600.;
const double MaxYDifLc  =  400.;

// TrackChiSqr
const double MaxChisqrSks    = 30.;
const double MaxChisqrBdcIn  = 30.;
const double MaxChisqrBdcOut = 30.;
const double MaxChisqrSdcIn  = 30.;
const double MaxChisqrSdcOut = 30.;


VEvent::VEvent()
{
}

VEvent::~VEvent()
{
}

class EventSksHodoscoepSigmaCut : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventSksHodoscoepSigmaCut();
  ~EventSksHodoscoepSigmaCut();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventSksHodoscoepSigmaCut::EventSksHodoscoepSigmaCut()
  : VEvent(), 
    rawData(0), 
    DCAna(new DCAnalyzer()), 
    hodoAna(new HodoAnalyzer())
{
}

EventSksHodoscoepSigmaCut::~EventSksHodoscoepSigmaCut()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 10
#endif

//For Tree
struct Event{
  int trigtype;
  int trigflag[NumOfMisc];

  double btof[MaxHits];

  int nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int ntIn;
  int nhIn[MaxHits]; 
  double chisqrIn[MaxHits];
  double u0In[MaxHits];
  double v0In[MaxHits];
  double x0In[MaxHits];
  double y0In[MaxHits];

  int nhx1;
  int x1wire[MaxHits];
  int nhx2;
  int x2wire[MaxHits];

  int ntOut;
  int nhOut[MaxHits]; 
  double chisqrOut[MaxHits];
  double u0Out[MaxHits];
  double v0Out[MaxHits];
  double x0Out[MaxHits];
  double y0Out[MaxHits];

  int ntSks;
  int nhSks[MaxHits]; 
  double chisqrSks[MaxHits];
  double path[MaxHits];
  double p[MaxHits];
  double m2[MaxHits];

  double xtgt[MaxHits];  
  double ytgt[MaxHits];  
  double utgt[MaxHits];  
  double vtgt[MaxHits]; 
  double theta[MaxHits]; 
  double phi[MaxHits]; 

  //Hodoscope
  int nhTof;
  int TofSeg[MaxHits];
  double tTof[MaxHits];
  double dtTof[MaxHits];
  double deTof[MaxHits];

  int nhLc;
  int LcSeg[MaxHits];
  double tLc[MaxHits];
  double dtLc[MaxHits];
  double deLc[MaxHits];

  int nhAc1;
  int Ac1Seg[NumOfSegAC];
  double tAc1[NumOfSegAC];
  double deAc1[NumOfSegAC];
  double npeAc1;

  int nhAc2;
  int Ac2Seg[NumOfSegAC];
  double tAc2[NumOfSegAC];
  double deAc2[NumOfSegAC];
  double npeAc2;
};
static Event event;

bool EventSksHodoscoepSigmaCut::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventSksHodoscoepSigmaCut::ProcessingNormal()
{
  static const std::string funcname = "[EventSksHodoscoepSigmaCut::ProcessingNormal]";

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
      //Bema x TOF
      //if( seg==3 && T>0 ){ trig_type = seg; event.trigtype = seg;}

      event.trigflag[seg] = T;
    }
  }

  //------------------------Cut
  if (trig_type != PIK_TRIG )  return true;
  //if (trig_type != BEAM_TRIG )  return true;
  //if (trig_type != PIPI_TRIG )  return true;
  //if (trig_type != BTOF_TRIG )  return true;
  HF1( 1, 1. );

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
  if( ncBh2!=1 ) return true;
  //if( ncBh2==0 ) return true;

  BH2Cluster *clBH2Time0=hodoAna->GetClusterBH2(0);
  double time0=clBH2Time0->CTime0();
  {
    int ncOk=0;
    double mint=clBH2Time0->CMeanTime();
    for( int i=0; i<ncBh2; ++i ){
      BH2Cluster *cl=hodoAna->GetClusterBH2(i);
      double t=cl->CMeanTime();
      double dEbh2 = cl->DeltaE();
      //------------------------Cut
      if( dEbh2<MinDeltaEBH2 || MaxDeltaEBH2<dEbh2 ) continue;
       event.Bh2Seg[i] = cl->MeanSeg()+1;
       event.tBh2[i]   = t;
       event.deBh2[i]  = cl->DeltaE();
       ++ncOk;
      if( fabs(t)<fabs(mint) ){
	clBH2Time0 = cl;
	mint=t; time0=clBH2Time0->CTime0();
      }
    }
    if( ncOk<1 ) return true;
  }

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData); 
  int ncBh1=hodoAna->GetNClustersBH1();
  event.nhBh1=ncBh1;
  {
    int ncOk=0;
    for( int i=0; i<ncBh1; ++i ){
      HodoCluster *cl=hodoAna->GetClusterBH1(i);
      double btof= (cl->CMeanTime())-time0;
      double dEbh1 = cl->DeltaE();
      //------------------------Cut
      if( dEbh1<MinDeltaEBH1 || MaxDeltaEBH1<dEbh1 ) continue;
      event.Bh1Seg[i]=cl->MeanSeg()+1;
      event.tBh1[i]=cl->CMeanTime();
      event.deBh1[i]=cl->DeltaE();
      event.btof[i] = btof;
      //------------------------Cut
      if( MinBeamToF<btof && btof<MaxBeamToF ){
	++ncOk;
      }
      else{
	cl->GoodForAnalysis( false );
      }
    }
    //------------------------Cut
    if( ncOk<1 ) return true;
  }

  HF1( 1, 2. );

  //////////////pre Event Selection by TOF
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();
  event.nhTof=ncTof;
  {
    int ncOk=0;
    for( int i=0; i<ncTof; ++i ){
      HodoCluster *cl=hodoAna->GetClusterTOF(i);
      int cs=cl->ClusterSize();
      double ms=cl->MeanSeg()+1, cmt=cl->CMeanTime(),
	de=cl->DeltaE(), dt=cl->TimeDif(),
	x=(ms-16.5)*70., y=dt*50.;
      double t=cmt-time0;	

      event.TofSeg[i]=ms;
      event.tTof[i]=t;
      event.dtTof[i]=dt;
      event.deTof[i]=de;

      //------------------------Cut
      if( MinDeltaETof<de && de<MaxDeltaETof &&
	  MinTimeTof  <t  && t<MaxTimeTof && (cl->MeanSeg()+1)>12 ){
        ++ncOk;
	HF1( TOFHid+31, double(cs) );
	HF1( TOFHid+32, ms-0.5 );
	HF1( TOFHid+33, cmt ); HF1( TOFHid+34, de ); HF1( TOFHid+35, dt );
	HF1( TOFHid+36, x ); HF1( TOFHid+37, y ); HF2( TOFHid+38, x, y );
	HF1( TOFHid+39, t );
      }
      else{
	cl->GoodForAnalysis( false );
      }
    }
    if( ncOk<1 ){
      //std::cout<<"#D Good TofCluster was not found. [return]"<<std::endl;
      return true;
    }
    HF1( TOFHid+30, double(ncTof) );
  }

  HF1( 1, 3. );

  //////////////pre Event Selection by LC
  hodoAna->DecodeLCHits(rawData);
  int ncLc=hodoAna->GetNClustersLC();
  event.nhLc=ncLc;
  {
    int ncOk=0;
    for( int i=0; i<ncLc; ++i ){
      HodoCluster *cl=hodoAna->GetClusterLC(i);
      int cs=cl->ClusterSize();
      double ms=cl->MeanSeg()+1, cmt=cl->CMeanTime(),
	de=cl->DeltaE(), dt=cl->TimeDif(),
	x=(ms-14.5)*100., y=dt*50.;
      double t=cl->CMeanTime()-time0;

      event.LcSeg[i]=ms;
      event.tLc[i]=t;
      event.dtLc[i]=dt;
      event.deLc[i]=de;

      //------------------------Cut
      if( MinDeltaELc<de && de<MaxDeltaELc &&
	  MinTimeLc  <t  && t<MaxTimeLc && (cl->MeanSeg()+1)>10 ){
        ++ncOk;
	HF1( LCHid+31, double(cs) );
	HF1( LCHid+32, ms-0.5 );
	HF1( LCHid+33, cmt ); HF1( LCHid+34, de ); HF1( LCHid+35, dt );
	HF1( LCHid+36, x ); HF1( LCHid+37, y ); HF2( LCHid+38, x, y );
	HF1( LCHid+39, t );
      }
      else{
        cl->GoodForAnalysis( false );
      }
    }
    if( ncOk<1 ){
      //std::cout<<"#D Good LcCluster was not found. [return]"<<std::endl;
      return true;
    }
    HF1( LCHid+30, double(ncLc) );
  }

  HF1( 1, 4. );
  
  //////////////pre Event Selection by TOF-LC
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
	//------------------------Cut
	if( MinToF2LSeg<=segTof-segLc && segTof-segLc<=MaxToF2LSeg
	    && MinToFT2L<tlc-ttof && tlc-ttof<MaxToFT2L ){
          ++nc;
	  HF2( 111, segTof-0.5, segLc-0.5 );
	  HF1( 113, segTof-segLc );
	  HF1( 115, tlc-ttof );
        }
        else{
	  clTof->GoodForAnalysis( false );
	  clLc->GoodForAnalysis( false );
        }
      }
    }
    if( nc<1 ){
      //std::cout<<"#D Good TOF-LC-Combination was not found. [return]"<<std::endl;
      return true;
    }
  }
  HF1( 1, 5. );
  
  //tree->Fill();return true;
  
  DCAna->DecodeRawHits( rawData );
  
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  int IdLc  = DCGeomMan::GetInstance().GetLcId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 

  double multi_SdcIn=0.;
  double multi_SdcOut=0.;

  //////////////SDC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn += double(nhIn);
    }
  }
  //------------------------Cut
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
   return true;

  //////////////SDC3&4 number of hit in one layer not 0
  //////////////SDC3-X1 vs SDC4-X2 chekc
  {
    int nhX1=0, nhX2=0;
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	int wire=hit->GetWire();
	if( layer == 2 ){ 
	  event.x1wire[nhX1] = wire;
	  nhX1++;
	}
	if( layer == 11 ){ 
	  event.x2wire[nhX2] = wire;
	  nhX2++;
	}
      }
    }
    event.nhx1=nhX1;
    event.nhx2=nhX2;
  }
  //------------------------Cut
  if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    return true;

  HF1( 1, 6. );

  //tree->Fill();return true;

  //////////////SdcIn tracking
  //std::cout<<"#D Begin TrackSearchSdcIn()"<<std::endl;
  DCAna->TrackSearchSdcIn();
  //std::cout<<"#D End TrackSearchSdcIn()"<<std::endl;

  int ntIn=DCAna->GetNtracksSdcIn();
  event.ntIn=ntIn;
  {
    int ntOk=0;
    for( int it=0; it<ntIn; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();
 
      HF1( 11, double(nh) );
      HF1( 12, chisqr ); 
      HF1( 14, x0 ); HF1( 15, y0 );
      HF1( 16, u0 ); HF1( 17, v0 );
      HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
      HF2( 20, x0, y0 );
      
      event.nhIn[it] = nh;
      event.chisqrIn[it] = chisqr;
      event.x0In[it] = x0;
      event.y0In[it] = y0;
      event.u0In[it] = u0;
      event.v0In[it] = v0;

      for( int ih=0; ih<nh; ++ih ){
	DCLTrackHit *hit=tp->GetHit(ih);
	int layerId=hit->GetLayer();
	HF1( 13, layerId );
      }
      //------------------------Cut
      if( chisqr<MaxChisqrSdcIn ) ++ntOk;
      else
        tp->GoodForTracking( false );
    }
    HF1( 10, double(ntOk) );
    if( ntOk<1 ) return true;
  }

  HF1( 1, 7. );

  //////////////SdcIn vs Tof&LC cut
  {
    int ntOk=0;
    for( int it=0; it<ntIn; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
      if(!tp) continue;
    
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX0(), y0=tp->GetY0();
      
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
      if( condTof && condLc ){
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
      else {
	tp->GoodForTracking( false );
      }
    }
    if( ntOk<1 ) return true;
  }

  HF1( 1, 8. );

  //////////////SdcOut tracking
  //std::cout<<"#D Begin TrackSearchSdcOut()"<<std::endl;
  DCAna->TrackSearchSdcOut();
  //std::cout<<"#D End TrackSearchSdcOut()"<<std::endl;
  
  int ntOut=DCAna->GetNtracksSdcOut();
  event.ntOut=ntOut;
  {
    int ntOk=0;
    for( int it=0; it<ntOut; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double x0=tp->GetX(zTof), y0=tp->GetY(zTof);
      double xlc=tp->GetX(zLc), ylc=tp->GetY(zLc);

      HF1( 31, double(nh) );
      HF1( 32, chisqr ); 
      HF1( 34, x0 ); HF1( 35, y0 );
      HF1( 36, u0 ); HF1( 37, v0 );
      HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
      HF2( 40, x0, y0 );
      
      event.nhOut[it] = nh;
      event.chisqrOut[it] = chisqr;
      event.x0Out[it] = x0;
      event.y0Out[it] = y0;
      event.u0Out[it] = u0;
      event.v0Out[it] = v0;
   
      bool condTof=false, condLc=false;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
        if( !clTof || !clTof->GoodForAnalysis() ) continue;
        double x=(clTof->MeanSeg()-15.5)*70.;
        double y=(clTof->TimeDif())*50.;
	HF1( TOFHid+41, x0-x ); HF1( TOFHid+42, y0-y );
	//------------------------Cut
  	if( MinXDifTof<x0-x && x0-x<MaxXDifTof &&
 	    MinYDifTof<y0-y && y0-y<MaxYDifTof )
          condTof=true;
      }
      for( int j=0; j<ncLc; ++j ){
        HodoCluster *clLc=hodoAna->GetClusterLC(j);
        if( !clLc || !clLc->GoodForAnalysis() ) continue;
        double x=(clLc->MeanSeg()-13.5)*100.;
        double y=(clLc->TimeDif())*50.;
	HF1( LCHid+41, xlc-x ); HF1( LCHid+42, ylc-y );
	//------------------------Cut
  	if( MinXDifLc<xlc-x && xlc-x<MaxXDifLc &&
	    MinYDifLc<ylc-y && ylc-y<MaxYDifLc )
          condLc=true;
      }
      if( condTof && condLc ){
	HF1( 31, double(nh) ); HF1( 32, chisqr ); 
	HF1( 34, x0 ); HF1( 35, y0 );
	HF1( 36, u0 ); HF1( 37, v0 );
	HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
	for( int ih=0; ih<nh; ++ih ){
	  DCLTrackHit *hit=tp->GetHit(ih);
	  int layerId=hit->GetLayer();
	  HF1( 33, layerId );
	}
        if( chisqr<MaxChisqrSdcOut ) ++ntOk;
	else 
	  tp->GoodForTracking( false );
      }
      else{
        tp->GoodForTracking( false );
      }
    }
    HF1( 30, double(ntOk) );
    if( ntOk<1 ) return true;
  }

  HF1( 1, 9. );

  //tree->Fill(); return true;

  //////////////SKS tracking
  //std::cout<<"#D Begin TrackSearchSks()"<<std::endl;
  DCAna->TrackSearchSks();
  //std::cout<<"#D End TrackSearchSks()"<<std::endl;

  int ntSks=DCAna->GetNTracksSks();
  event.ntSks = ntSks; 
  HF1( 50, double(ntSks) );
  {
    int ntOk=0;
    for( int i=0; i<ntSks; ++i ){
      SksTrack *tp=DCAna->GetSksTrack(i);
      if(!tp) continue;
      int nh=tp->GetNHits();
      double chisqr=tp->chisqr();
      ThreeVector Ppos=tp->PrimaryPosition();
      ThreeVector Pmom=tp->PrimaryMomentum();
      double pathL=tp->PathLengthToTOF();
      double xt=Ppos.x(), yt=Ppos.y();
      double p=Pmom.mag();
      double ut=Pmom.x()/p, vt=Pmom.y()/p;
      double cost = 1./sqrt(1.+ut*ut+vt*vt);
      double theta = acos(cost)*Rad2Deg;
      double phi   = atan2( ut, vt );

      HF1( 51, double(nh) );
      HF1( 52, chisqr );
      HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
      HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
      HF1( 61, p ); HF1( 62, pathL );

      event.nhSks[i] = nh;
      event.chisqrSks[i] = chisqr;
      event.path[i] = pathL;
      event.p[i] = p;
      
      event.xtgt[i] = xt;
      event.ytgt[i] = yt;
      event.utgt[i] = ut;
      event.vtgt[i] = vt; 
      
      event.theta[i] = theta;

      std::cout<<"******************************"<<std::endl;
      double m2;
      for( int j=0; j<ncTof; ++j ){
	HodoCluster *clTof=hodoAna->GetClusterTOF(j);
	if( !clTof ) continue;
	double time = clTof->CMeanTime()-time0+offset;
	
	if( time>0 ){
	  HF1( 63, m2 );
	  m2=MassSquare( p, pathL, time );
// 	  std::cout<<"Mom= "<< p <<std::endl;
// 	  std::cout<<"Path= "<< pathL <<std::endl;
// 	  std::cout<<"Time= "<< time <<std::endl;
 	  std::cout<<"m2= "<< m2 <<std::endl;
	}
      }      
      event.m2[i] = m2;

      for( int j=0; j<nh; ++j ){
	TrackHit *hit=tp->GetHit(j);
	if(!hit) continue;
	int layerId=hit->GetLayer();
	HF1(53, layerId );
      }
      //------------------------Cut
      if( chisqr<MaxChisqrSks ) ++ntOk;
      else continue;
      
      DCLocalTrack *trIn =tp->GetLocalTrackIn();
      DCLocalTrack *trOut=tp->GetLocalTrackOut();
      if(trIn){
	int nhIn=trIn->GetNHit();
	double x0in=trIn->GetX0(), y0in=trIn->GetY0();
	double u0in=trIn->GetU0(), v0in=trIn->GetV0();
	double chiin=trIn->GetChiSquare();
	HF1( 21, double(nhIn) ); HF1( 22, chiin );
	HF1( 24, x0in ); HF1( 25, y0in ); HF1( 26, u0in ); HF1( 27, v0in );
	HF2( 28, x0in, u0in ); HF2( 29, y0in, v0in );
	for( int jin=0; jin<nhIn; ++jin ){
	  int layer=trIn->GetHit(jin)->GetLayer();
	  HF1( 23, layer );
	}
      }
      if(trOut){
	int nhOut=trOut->GetNHit();
	double x0out=trOut->GetX(zTof), y0out=trOut->GetY(zTof);
	double u0out=trOut->GetU0(), v0out=trOut->GetV0();
	double chiout=trOut->GetChiSquare();
	HF1( 41, double(nhOut) ); HF1( 42, chiout );
	HF1( 44, x0out ); HF1( 45, y0out ); HF1( 46, u0out ); HF1( 47, v0out );
	HF2( 48, x0out, u0out ); HF2( 49, y0out, v0out );
	for( int jout=0; jout<nhOut; ++jout ){
	  int layer=trOut->GetHit(jout)->GetLayer();
	  HF1( 43, layer );
	}
	//HF2( 64, x0out, p );
      }
      //HF1( 50, double(ntOk) );
      if( ntOk<1 ) return true;
    }
  }

  HF1( 1, 10. );

  //------------------------Cut
  if( ntSks!=1 ) return true;

  HF1( 1, 11. );

  SksTrack *track=DCAna->GetSksTrack(0);
  double plTof=track->PathLengthToTOF();
  double plLc=track->PathLengthTotal();
  double p=track->PrimaryMomentum().mag();
  double tpiTof=plTof*sqrt(p*p+PionMass*PionMass)/p/LightVel;
  double tpiLc=plLc*sqrt(p*p+PionMass*PionMass)/p/LightVel;
  double tkTof=plTof*sqrt(p*p+KaonMass*KaonMass)/p/LightVel;
  double tkLc=plLc*sqrt(p*p+KaonMass*KaonMass)/p/LightVel;
  double tpTof=plTof*sqrt(p*p+ProtonMass*ProtonMass)/p/LightVel;
  double tpLc=plLc*sqrt(p*p+ProtonMass*ProtonMass)/p/LightVel;

  HF1( 71, plTof ); HF1( 72, plLc ); HF1( 73, p );
  HF1( 74, tpiTof ); HF1( 75, tpiLc );
  HF1( 76, tkTof ); HF1( 77, tkLc );
  HF1( 78, tpTof ); HF1( 79, tpLc );

  // TOF
  {
    int nh=hodoAna->GetNHitsTOF();
    HF1( TOFHid+50, double(nh) );
     for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitTOF(i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( TOFHid+51, seg-0.5 );
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double cmt=hit->CMeanTime(), t= cmt-time0+offset;//cmt-time0;
      double de=hit->DeltaE();
      double m2=MassSquare( p, plTof, t );
      double beta=plTof/t/LightVel;
      HF1( TOFHid+100*seg+31, cmt );
      HF1( TOFHid+100*seg+32, ctd-ctu );
      HF1( TOFHid+100*seg+33, t );
      HF1( TOFHid+100*seg+34, t-tpiTof );
      HF1( TOFHid+100*seg+35, t-tkTof );
      HF1( TOFHid+100*seg+36, t-tpTof );
      HF1( TOFHid+100*seg+37, m2 );
      HF1( TOFHid+100*seg+38, beta );
      HF1( TOFHid+52, cmt ); HF1( TOFHid+53, de );
      HF1( TOFHid+54, ctd-ctu ); HF1( TOFHid+55, t );
      HF1( TOFHid+56, t-tpiTof ); HF1( TOFHid+57, t-tkTof );
      HF1( TOFHid+58, t-tpTof ); HF1( TOFHid+59, m2 );
      HF1( TOFHid+60, beta );
    }
  }

  // LC
  {
    int nh=hodoAna->GetNHitsLC();
    HF1( LCHid+50, double(nh) );
    int nh2=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitLC(i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( LCHid+51, seg-0.5 );
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double cmt=hit->CMeanTime(), t= cmt-time0+offset;//cmt-time0;
      double de=hit->DeltaE();
      double m2=MassSquare( p, plLc, t );
      double beta=plLc/t/LightVel;
      HF1( LCHid+100*seg+31, cmt );
      HF1( LCHid+100*seg+32, ctd-ctu );
      HF1( LCHid+100*seg+33, t );
      HF1( LCHid+100*seg+34, t-tpiLc );
      HF1( LCHid+100*seg+35, t-tkLc );
      HF1( LCHid+100*seg+36, t-tpLc );
      HF1( LCHid+100*seg+37, m2 );
      HF1( LCHid+100*seg+38, beta );
      HF1( LCHid+52, cmt ); HF1( LCHid+53, de );
      HF1( LCHid+54, ctd-ctu ); HF1( LCHid+55, t );
      HF1( LCHid+56, t-tpiLc ); HF1( LCHid+57, t-tkLc );
      HF1( LCHid+58, t-tpLc ); HF1( LCHid+59, m2 );
      HF1( LCHid+60, beta );
    }
  }

  // AC
  hodoAna->DecodeACHits(rawData);
  for( int layer=1; layer<=2; ++layer ){
    int nh=hodoAna->GetNHitsAC(layer);
    HF1( ACHid+10+50*(layer-1), double(nh) );
    if( layer==1 ) event.nhAc1 = nh;
    if( layer==2 ) event.nhAc2 = nh;
    int nh2=0;
    double npeac1, npeac2;
    for( int i=0; i<nh; ++i ){
      Hodo1Hit *hit=hodoAna->GetHitAC(layer,i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( ACHid+11+50*(layer-1), seg-0.5 );
      double a=hit->GetA(), t=hit->GetT(), ct=hit->GetCT();
      HF1( ACHid+12+50*(layer-1), a );
      HF1( ACHid+13+50*(layer-1), t );

      if( layer==1 ){
	event.Ac1Seg[i] = seg;
	event.deAc1[i] = a;
	event.tAc1[i] = t;
	if( t >0 ) npeac1 += a;
      }
      if( layer==2 ){
	event.Ac2Seg[i] = seg;
	event.deAc2[i] = a;
	event.tAc2[i] = t;
	if( t >0 ) npeac2 += a;
      }
      if(a>0.1)	++nh2;
    }
    event.npeAc1 = npeac1;
    event.npeAc2 = npeac2;

    HF1( ACHid+20+50*(layer-1), double(nh2) );
  }

  // TOF-LC-AC
  {
    int nhTof=hodoAna->GetNHitsTOF();
    int nhLc =hodoAna->GetNHitsLC();
    int nhAc1=hodoAna->GetNHitsAC(1);
    int nhAc2=hodoAna->GetNHitsAC(2);
    for( int iTof=0; iTof<nhTof; ++iTof ){
      Hodo2Hit *hitTof=hodoAna->GetHitTOF(iTof);
      if(!hitTof) continue;
      int segTof=hitTof->SegmentId()+1; 
      double cmtTof=hitTof->CMeanTime();
      double deTof=hitTof->DeltaE();
      double m2=MassSquare( p, plTof, cmtTof-time0+offset );
      HF2( 301, m2, deTof );
      for( int iLc=0; iLc<nhLc; ++iLc ){
        Hodo2Hit *hitLc =hodoAna->GetHitLC(iLc);
        if(!hitLc) continue;
        int segLc=hitLc->SegmentId()+1;
        double cmtLc=hitLc->CMeanTime();
        double deLc=hitLc->DeltaE();
        HF2( 201, segTof-0.5, segLc-0.5 );
        HF1( 203, cmtLc-cmtTof );
	HF2( 302, m2, deLc );
	HF2( 303, m2, cmtLc-cmtTof );
      }
      for( int iAc=0; iAc<nhAc1; ++iAc ){
	Hodo1Hit *hit=hodoAna->GetHitAC(1,iAc);
	if(!hit) continue;
	double a=hit->GetA();
	HF2( 304, m2, a);
      }
      for( int iAc=0; iAc<nhAc2; ++iAc ){
	Hodo1Hit *hit=hodoAna->GetHitAC(2,iAc);
	if(!hit) continue;
	double a=hit->GetA();
	HF2( 305, m2, a);
      }
    }
  }
  
  HF1( 1, 12. );
  tree->Fill();
  //std::cout<<"******************************"<<std::endl;

  return true;
}

void EventSksHodoscoepSigmaCut::InitializeEvent( void )
{

  event.trigtype = -1;  
  event.nhBh2  = -1;
  event.nhBh1  = -1;
  event.ntIn     = -1;
  event.ntOut    = -1;
  event.ntSks    = -1;
  event.nhx1     = -1;
  event.nhx2     = -1;

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it] = -9999.0;
    event.deBh2[it] = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it] = -9999.0;
    event.deBh1[it] = -9999.0;
    event.btof[it]  = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.x1wire[it] = -1;
    event.x2wire[it] = -1;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhIn[it]     = -1;
    event.chisqrIn[it] = -1.0;
    event.x0In[it] = -9999.0;
    event.y0In[it] = -9999.0;
    event.u0In[it] = -9999.0;
    event.v0In[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhOut[it]     = -1; 
    event.chisqrOut[it] = -1.0;
    event.x0Out[it] = -9999.0;
    event.y0Out[it] = -9999.0;
    event.u0Out[it] = -9999.0;
    event.v0Out[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSks[it]     = -1; 
    event.chisqrSks[it] = -1.0;
    event.path[it]  = -9999.0;
    event.p[it]     = -9999.0;
    event.m2[it]    = -9999.0;

    event.xtgt[it]  = -9999.0;
    event.ytgt[it]  = -9999.0;
    event.utgt[it]  = -9999.0;
    event.vtgt[it]  = -9999.0;
    event.theta[it] = -9999.0;
    event.phi[it]   = -9999.0;
  }

  //Hodoscope
  event.nhTof  = -1;
  event.nhLc   = -1;
  event.nhAc1  = -1;
  event.nhAc2  = -1;
  event.npeAc1  = -9999.0;
  event.npeAc2  = -9999.0;

  for( int it=0; it<MaxHits; it++){
    event.TofSeg[it] = -1;
    event.tTof[it] = -9999.0;
    event.dtTof[it] = -9999.0;
    event.deTof[it] = -9999.0;

    event.LcSeg[it] = -1;
    event.tLc[it] = -9999.0;
    event.dtLc[it] = -9999.0;
    event.deLc[it] = -9999.0;

    event.Ac1Seg[it] = -1;
    event.tAc1[it] = -9999.0;
    event.deAc1[it] = -9999.0;

    event.Ac2Seg[it] = -1;
    event.tAc2[it] = -9999.0;
    event.deAc2[it] = -9999.0;
  }
}

bool EventSksHodoscoepSigmaCut::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSksHodoscoepSigmaCut;
}

bool ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );

  HB1( 10, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 11, "#Hits of Track SdcIn", 15, 0., 15. );
  HB1( 12, "Chisqr SdcIn", 500, 0., 50. ); 
  HB1( 13, "LayerId SdcIn", 15, 0., 15. );
  HB1( 14, "X0 SdcIn", 400, -100., 100. ); 
  HB1( 15, "Y0 SdcIn", 400, -100., 100. );
  HB1( 16, "U0 SdcIn", 200, -0.20, 0.20 );
  HB1( 17, "V0 SdcIn", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100., 100. );

  HB1( 21, "#Hits of Track SdcIn [SksTrack]", 15, 0., 15. );
  HB1( 22, "Chisqr SdcIn [SksTrack]", 500, 0., 50. ); 
  HB1( 23, "LayerId SdcIn [SksTrack]", 15, 0., 15. );
  HB1( 24, "X0 SdcIn [SksTrack]", 400, -100., 100. ); 
  HB1( 25, "Y0 SdcIn [SksTrack]", 400, -100., 100. );
  HB1( 26, "U0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
  HB1( 27, "V0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
  HB2( 28, "U0%X0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 29, "V0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
  //HB2( 30, "X0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -100., 100. );

  HB1( 30, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 31, "#Hits of Track SdcOut", 20, 0., 20. );
  HB1( 32, "Chisqr SdcOut", 500, 0., 50. ); 
  HB1( 33, "LayerId SdcOut", 20, 30., 50. );
  HB1( 34, "X0 SdcOut", 1400, -1200., 1200. ); 
  HB1( 35, "Y0 SdcOut", 1000, -500., 500. );
  HB1( 36, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 37, "V0 SdcOut",  200, -0.20, 0.20 );
  HB2( 38, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40 );
  HB2( 39, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20 );
  HB2( 40, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500., 500. );

  HB1( 41, "#Hits of Track SdcOut [SksTrack]", 20, 0., 20. );
  HB1( 42, "Chisqr SdcOut [SksTrack]", 500, 0., 50. ); 
  HB1( 43, "LayerId SdcOut [SksTrack]", 20, 30., 50. );
  HB1( 44, "X0 SdcOut [SksTrack]", 1400, -1200., 1200. ); 
  HB1( 45, "Y0 SdcOut [SksTrack]", 1000, -500., 500. );
  HB1( 46, "U0 SdcOut [SksTrack]",  700, -0.35, 0.35 );
  HB1( 47, "V0 SdcOut [SksTrack]",  200, -0.10, 0.10 );
  HB2( 48, "U0%X0 SdcOut [SksTrack]", 120, -600., 600., 100, -0.40, 0.40 );
  HB2( 49, "V0%Y0 SdcOut [SksTrack]", 100, -500., 500., 100, -0.10, 0.10 );
  //HB2( 50, "X0%Y0 SdcOut [SksTrack]", 100, -700., 700., 100, -500., 500. );

  HB1( 50, "#Tracks SKS", 10, 0., 10. );
  HB1( 51, "#Hits of SksTrack", 30, 0., 30. );
  HB1( 52, "Chisqr SksTrack", 500, 0., 100. );
  HB1( 53, "LayerId SksTrack", 50, 0., 50. );
  HB1( 54, "Xtgt SksTrack", 200, -100., 100. );
  HB1( 55, "Ytgt SksTrack", 200, -100., 100. );
  HB1( 56, "Utgt SksTrack", 300, -0.30, 0.30 );
  HB1( 57, "Vtgt SksTrack", 300, -0.20, 0.20 );
  HB2( 58, "U%Xtgt SksTrack", 100, -100., 100., 100, -0.25, 0.25 );
  HB2( 59, "V%Ytgt SksTrack", 100, -100., 100., 100, -0.10, 0.10 );
  HB2( 60, "Y%Xtgt SksTrack", 100, -100., 100., 100, -100., 100. );
  HB1( 61, "P SksTrack", 400, 0.50, 0.90 );
  HB1( 62, "PathLength SksTrack", 600, 3000., 6000. );
  HB1( 63, "MassSqr", 600, -0.2, 1.2 );
  //   HB1( 62, "PathLength SksTrack", 400, 3000., 7000. );
  //   HB2( 63, "P%PathLength SksTrack", 100, 3000., 7000., 100, 0.60, 0.90 );
  //   HB2( 64, "P%X(zTof) SksTrack", 200, -1200., 1200., 100, 0.60, 0.90 );

  HB1( 71, "PathLength to TOF", 400, 3000., 7000. );
  HB1( 72, "PathLength to LC", 400, 3000., 7000. );
  HB1( 73, "P", 300, 0.60, 0.90 );
  HB1( 74, "Tpi Tof", 500, -5., 45. );
  HB1( 75, "Tpi Lc", 500, -5., 45. );
  HB1( 76, "Tk Tof", 500, -5., 45. );
  HB1( 77, "Tk Lc", 500, -5., 45. );
  HB1( 78, "Tp Tof", 500, -5., 45. );
  HB1( 79, "Tp Lc", 500, -5., 45. );

  // BH2
  HB1( BH2Hid+10, "#Hits BH2[Hodo]", NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+50, "CTime0 BH2", 200, -10., 10.);

  // Tof
  HB1( TOFHid+30, "#Clusters Tof", NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+31, "ClusterSize Tof", 5, 0., 5. ); 
  HB1( TOFHid+32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., double(NumOfSegTOF) );
  HB1( TOFHid+33, "CMeamTime Cluster Tof", 500, -5., 45. ); 
  HB1( TOFHid+34, "DeltaE Cluster Tof", 200, -0.5, 4.5 );
  HB1( TOFHid+35, "TimeDif Cluster Tof", 200, -10., 10. );
  HB1( TOFHid+36, "PosX Cluster Tof", 2*NumOfSegTOF, -double(NumOfSegTOF)*35., double(NumOfSegTOF)*35. );
  HB1( TOFHid+37, "PosY Cluster Tof", 200, -600., 600. );
  HB2( TOFHid+38, "PosY%PosX Cluster Tof", 2*NumOfSegTOF, -double(NumOfSegTOF)*35., double(NumOfSegTOF)*35., 200, -600., 600. );
  HB1( TOFHid+39, "CMeamTime-Time0 Cluster Tof", 500, -5., 45. ); 

  HB1( TOFHid+41, "X(zTof)-PosXTof", 500, -500., 500. );
  HB1( TOFHid+42, "Y(zTof)-PosYTof", 500, -500., 500. );

  HB1( TOFHid+50, "#Hits Tof[SksHodo]", NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+51, "Hitpat Tof[SksHodo]", NumOfSegTOF, 0., double(NumOfSegTOF) );
  HB1( TOFHid+52, "CMeanTime Tof[SksHodo]", 500, -5., 45. );
  HB1( TOFHid+53, "dE Tof[SksHodo]", 200, -0.5, 4.5 );
  HB1( TOFHid+54, "Tdown-Tup Tof[SksHodo]", 200, -10.0, 10.0 );
  HB1( TOFHid+55, "TofTime[SksHodo]", 500, -5., 45. );
  HB1( TOFHid+56, "TofTime-Tpi", 500, -25., 25. );
  HB1( TOFHid+57, "TofTime-Tk", 500, -25., 25. );
  HB1( TOFHid+58, "TofTime-Tp", 500, -25., 25. );
  HB1( TOFHid+59, "MassSquare Tof", 600, -0.1, 1.1 );
  HB1( TOFHid+60, "beta Tof", 500, 0., 1.5 );

  for( int i=1; i<=NumOfSegTOF; ++i ){
    std::stringstream title1,title2;
    title1 << "Tof-" << i << " CMeanTime[SksHodo]";
    HB1( TOFHid+100*i+31, title1.str().c_str(), 500, -5., 45. );
    title2 << "TOF-" << i << " Tdown-Tup[SksHodo]";
    HB1( TOFHid+100*i+32, title2.str().c_str(), 200, -10.0, 10.0 );
    std::stringstream title3,title4,title5,title6;
    title3 << "Tof-" << i << " TofTime[SksHodo]";
    HB1( TOFHid+100*i+33, title3.str().c_str(), 500, -5., 45. );
    title4 << "Tof-" << i << " TofTime-Tpi";
    HB1( TOFHid+100*i+34, title4.str().c_str(), 500, -25., 25. );
    title5 << "Tof-" << i << " TofTime-Tk";
    HB1( TOFHid+100*i+35, title5.str().c_str(), 500, -25., 25. );
    title6 << "Tof-" << i << " TofTime-Tp";
    HB1( TOFHid+100*i+36, title6.str().c_str(), 500, -25., 25. );
    std::stringstream title7,title8;
    title7 << "Tof-" << i << " MassSquare";
    HB1( TOFHid+100*i+37, title7.str().c_str(), 600, -0.1, 1.1 );
    title8 << "Tof-" << i << " beta";
    HB1( TOFHid+100*i+38, title8.str().c_str(), 500, 0., 1.5 );
  }

  // Lc
  HB1( LCHid+30, "#Clusters Lc", NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+31, "ClusterSize Lc", 5, 0., 5. ); 
  HB1( LCHid+32, "HitPat Cluster Lc", 2*NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( LCHid+33, "CMeamTime Cluster Lc", 500, -5., 45. ); 
  HB1( LCHid+34, "DeltaE Cluster Lc", 200, -0.5, 4.5 );
  HB1( LCHid+35, "TimeDif Cluster Lc", 200, -15., 15. );
  HB1( LCHid+36, "PosX Cluster Lc", 2*NumOfSegLC, -double(NumOfSegLC)*50., double(NumOfSegLC)*50. );
  HB1( LCHid+37, "PosY Cluster Lc", 200, -800., 800. );
  HB2( LCHid+38, "PosY%PosX Cluster Lc", 2*NumOfSegLC, -double(NumOfSegLC)*50., double(NumOfSegLC)*50., 200, -800., 800. );
  HB1( LCHid+39, "CMeamTime-Time0 Cluster Lc", 500, -5., 45. ); 

  HB1( LCHid+41, "X(zLc)-PosXLc", 500, -500., 500. );
  HB1( LCHid+42, "Y(zLc)-PosYLc", 500, -500., 500. );

  HB1( LCHid+50, "#Hits Lc[SksHodo]", NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+51, "Hitpat Lc[SksHodo]", NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( LCHid+52, "CMeanTime Lc[SksHodo]", 500, -5., 45. );
  HB1( LCHid+53, "dE Lc[SksHodo]", 200, -0.5, 4.5 );
  HB1( LCHid+54, "Tdown-Tup Lc[SksHodo]", 200, -10.0, 10.0 );
  HB1( LCHid+55, "LcTime[SksHodo]", 500, -5., 45. );
  HB1( LCHid+56, "LcTime-Tpi", 500, -25., 25. );
  HB1( LCHid+57, "LcTime-Tk", 500, -25., 25. );
  HB1( LCHid+58, "LcTime-Tp", 500, -25., 25. );
  HB1( LCHid+59, "MassSquare Lc", 600, -0.1, 1.1 );
  HB1( LCHid+60, "beta Lc", 500, 0., 1.5 );

  for( int i=1; i<=NumOfSegLC; ++i ){
    std::stringstream title1,title2;
    title1 << "Lc-" << i << " CMeanTime[SksHodo]";
    HB1( LCHid+100*i+31, title1.str().c_str(), 500, -5., 45. );
    title2 << "Lc-" << i << " Tdown-Tup[SksHodo]";
    HB1( LCHid+100*i+32, title2.str().c_str(), 200, -10.0, 10.0 );
    std::stringstream title3,title4,title5,title6;
    title3 << "Lc-" << i << " LcTime[SksHodo]";
    HB1( LCHid+100*i+33, title3.str().c_str(), 500, -5., 45. );
    title4 << "Lc-" << i << " LcTime-Tpi";
    HB1( LCHid+100*i+34, title4.str().c_str(), 500, -25., 25. );
    title5 << "Lc-" << i << " LcTime-Tk";
    HB1( LCHid+100*i+35, title5.str().c_str(), 500, -25., 25. );
    title6 << "Lc-" << i << " LcTime-Tp";
    HB1( LCHid+100*i+36, title6.str().c_str(), 500, -25., 25. );
    std::stringstream title7,title8;
    title7 << "Lc-" << i << " MassSquare";
    HB1( LCHid+100*i+37, title7.str().c_str(), 600, -0.1, 1.1 );
    title8 << "Lc-" << i << " beta";
    HB1( LCHid+100*i+38, title8.str().c_str(), 500, 0., 1.5 );
  }

  // Ac
  HB1( ACHid+10, "#Hits Ac1[Hodo]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );
  HB1( ACHid+11, "Hitpat Ac1[Hodo]", NumOfSegAC,   0., double(NumOfSegAC)   );
  HB1( ACHid+12, "dE Ac1", 200, -0.5, 4.5 );
  HB1( ACHid+13, "T Ac1", 200,  0.0, 100 );
  HB1( ACHid+20, "#Hits Ac1[HodoGood]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );

  HB1( ACHid+60, "#Hits Ac2[Hodo]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );
  HB1( ACHid+61, "Hitpat Ac2[Hodo]", NumOfSegAC,   0., double(NumOfSegAC)   );
  HB1( ACHid+62, "dE Ac2", 200, -0.5, 4.5 );
  HB1( ACHid+63, "T Ac2", 200,  0.0, 100 );
  HB1( ACHid+70, "#Hits Ac2[HodoGood]",  NumOfSegAC+1, 0., double(NumOfSegAC+1) );

  // TOF-LC-AC
  HB2( 111, "SegLc%SegTof Cluster", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( 113, "SegTof-SegLc Cluster", 30, -15., 15. );
  HB1( 115, "LcTime-TofTime Cluster", 300, -5., 25. );

  HB2( 201, "SegLc%SegTof[SksHodo]", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( 203, "LcTime-TofTime[SksHodo]", 300, -5., 25. );

  HB2( 301, "TofDe%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );
  HB2( 302, "LcDe%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );
  HB2( 303, "LcTime-TofTime%MSqrTof", 120, -0.1, 1.1, 100, -5., 25 );
  HB2( 304, "Ac1De%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );
  HB2( 305, "Ac2De%MSqrTof", 120, -0.1, 1.1, 100, -0.5, 4.5 );

 ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");

  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

  tree->Branch("nhx1",    &event.nhx1,     "nhx1/I");
  tree->Branch("x1wire",   event.x1wire,   "x1wire[nhx1]/I");
  tree->Branch("nhx2",    &event.nhx2,     "nhx2/I");
  tree->Branch("x2wire",   event.x2wire,   "x2wire[nhx2]/I");

  tree->Branch("ntIn",    &event.ntIn,     "ntIn/I");
  tree->Branch("nhIn",     event.nhIn,     "nhIn[ntIn]/I");
  tree->Branch("chisqrIn", event.chisqrIn, "chisqrIn[ntIn]/D");
  tree->Branch("x0In",     event.x0In,     "x0In[ntIn]/D");
  tree->Branch("y0In",     event.y0In,     "y0In[ntIn]/D");
  tree->Branch("u0In",     event.u0In,     "u0In[ntIn]/D");
  tree->Branch("v0In",     event.v0In,     "v0In[ntIn]/D");

  tree->Branch("ntOut",   &event.ntOut,     "ntOut/I");
  tree->Branch("nhOut",    event.nhOut,     "nhOut[ntOut]/I");
  tree->Branch("chisqrOut",event.chisqrOut, "chisqrIn[ntOut]/D");
  tree->Branch("x0Out",    event.x0Out,     "x0Out[ntOut]/D");
  tree->Branch("y0Out",    event.y0Out,     "y0Out[ntOut]/D");
  tree->Branch("u0Out",    event.u0Out,     "u0Out[ntOut]/D");
  tree->Branch("v0Out",    event.v0Out,     "v0Out[ntOut]/D");

  tree->Branch("ntSks",      &event.ntSks,     "ntSks/I");
  tree->Branch("nhSks",       event.nhSks,     "nhSks[ntSks]/I");
  tree->Branch("chisqrSks",   event.chisqrSks, "chisqrSks[ntSks]/D");
  tree->Branch("path",        event.path,      "path[ntSks]/D");
  tree->Branch("p",           event.p,         "p[ntSks]/D");
  tree->Branch("m2",          event.m2,        "m2[ntSks]/D");

  tree->Branch("xtgt",    event.xtgt,   "xtgt[ntSks]/D");
  tree->Branch("ytgt",    event.ytgt,   "ytgt[ntSks]/D");
  tree->Branch("utgt",    event.utgt,   "utgt[ntSks]/D");
  tree->Branch("vtgt",    event.vtgt,   "vtgt[ntSks]/D");

  tree->Branch("theta",   event.theta,  "theta[ntSks]/D");
  tree->Branch("phi",     event.phi,    "phi[ntSks]/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/I");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhLc",   &event.nhLc,   "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/I");
  tree->Branch("tLc",     event.tLc,    "tLc[nhLc]/D");
  tree->Branch("dtLc",    event.dtLc,   "dtLc[nhLc]/D");
  tree->Branch("deLc",    event.deLc,   "deLc[nhLc]/D");

  tree->Branch("nhAc1",   &event.nhAc1,   "nhAc1/I");
  tree->Branch("Ac1Seg",   event.Ac1Seg,  "Ac1Seg[nhAc1]/I");
  tree->Branch("tAc1",     event.tAc1,    "tAc1[nhAc1]/D");
  tree->Branch("deAc1",    event.deAc1,   "deAc1[nhAc1]/D");
  tree->Branch("npeAc1",  &event.npeAc1,  "npeAc1/D");

  tree->Branch("nhAc2",   &event.nhAc2,   "nhAc2/I");
  tree->Branch("Ac2Seg",   event.Ac2Seg,  "Ac2Seg[nhAc2]/I");
  tree->Branch("tAc2",     event.tAc2,    "tAc2[nhAc2]/D");
  tree->Branch("deAc2",    event.deAc2,   "deAc2[nhAc2]/D");
  tree->Branch("npeAc2",  &event.npeAc2,  "npeAc2/D");

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

  return true;
}



