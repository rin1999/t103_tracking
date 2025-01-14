/*
  UserK18Tracking.cc
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

const double MaxMultiHitBcIn  =  5.0;
const double MaxMultiHitBcOut =  100.;

#define In  1
#define Out 1
#define K18 1

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventK18Tracking 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventK18Tracking();
  ~EventK18Tracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventK18Tracking::EventK18Tracking()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventK18Tracking::~EventK18Tracking()
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
  int trigtype;
  int trigflag[NumOfMisc];

  //Beam
 double btof[MaxHits];

  int nhBh2;
  double Bh2Seg[MaxHits];
  double tBh2[MaxHits];
  double deBh2[MaxHits];

  int nhBh1;
  double Bh1Seg[MaxHits];
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  //BcIn
  int ntIn;
  int nhIn[MaxHits]; 
  double chisqrIn[MaxHits];
  double u0In[MaxHits];
  double v0In[MaxHits];
  double x0In[MaxHits];
  double y0In[MaxHits];

  //BcOut
  int ntOut;
  int nhOut[MaxHits]; 
  double chisqrOut[MaxHits];
  double u0Out[MaxHits];
  double v0Out[MaxHits];
  double x0Out[MaxHits];
  double y0Out[MaxHits];

  //K1.8
  int ntK18;
  int nhK18[MaxHits]; 
  double chisqrK18[MaxHits];
  double p[MaxHits];
  double Delta[MaxHits];

  double xtgt[MaxHits];  
  double ytgt[MaxHits];  
  double utgt[MaxHits];  
  double vtgt[MaxHits]; 
  double theta[MaxHits];
  double phi[MaxHits]; 
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


bool EventK18Tracking::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventK18Tracking::ProcessingNormal()
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
      if( seg==1 && T>0 ){ trig_type = seg; event.trigtype = seg;}
      //(pi, K)
      //if( seg==2 && T>0 ){ trig_type = seg; event.trigtype = seg;}
      //(pi, pi)
      //if( seg==11 && T>0 ){ trig_type = seg; event.trigtype = seg;}

      event.trigflag[seg] = T;
    }
  }

  //------------------------Cut
  //if (trig_type != PIK_TRIG )  return true;
  if (trig_type != BEAM_TRIG )  return true;
  //if (trig_type != PIPI_TRIG )  return true;

  HF1( 1, 1. );

  //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
  //if( ncBh2!=1 ) return true;
  if( ncBh2==0 ) return true;

  HF1( 1, 2. );

  //////////////BH2 Analysis
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

  HF1( 1, 3. );

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
  HF1( 1, 4. );

  //tree->Fill(); return true;

  HF1( 1, 10. );

  DCAna->DecodeRawHits( rawData );  

  //////////////BC3&4 number of hit in one layer not 0
  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    return true;

  HF1( 1, 11. );

  //////////////BCOut tracking
#if Out
  DCAna->TrackSearchBcOut();

  int ntOut=DCAna->GetNtracksBcOut(); 
  event.ntOut=ntOut;
  int ntOutOk=0;
  HF1( 30, double(ntOut) );
  for( int it=0; it<ntOut; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX(0.), y0=tp->GetY(0.);

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

    if( chisqr < 20 ) ntOutOk++;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-100;
      HF1( 33, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+1, wire-0.5 );
      HF1( 100*layerId+2, dt );
      HF1( 100*layerId+3, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+4, pos );
      HF1( 100*layerId+5, res );
      HF2( 100*layerId+6, pos, res );
      HF2( 100*layerId+7, xcal, ycal);
    }
  }

//   for( int i1=0; i1<ntIn; ++i1 ){
//     DCLocalTrack *trIn=DCAna->GetTrackBcIn(i1);
//     double yin=trIn->GetY(500.), vin=trIn->GetV0();
//     for( int i2=0; i2<ntOut; ++i2 ){
//       DCLocalTrack *trOut=DCAna->GetTrackBcOut(i2);
//       double yout=trOut->GetY(3800.), vout=trOut->GetV0();
//       HF2( 20001, yin, yout ); HF2( 20002, vin, vout );
//       HF2( 20003, vin, yout ); HF2( 20004, vout, yin );
//     }
//   }
  //if( ntOut !=1 ) return true;
  //if( ntOut<1 ) return true;
  if( ntOutOk<1 ) return true;
#endif
  HF1( 1, 12. );

  HF1( 1, 15. );

  //////////////BC1&2 number of hit layer
  double multi_BcIn=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcInHC(layer);
      int nhIn=contIn.size();
      multi_BcIn += double(nhIn);
    }
  }
  if( multi_BcIn/double(NumOfLayersBcIn) > MaxMultiHitBcIn )
   return true;

  HF1( 1, 16. );

  //////////////BCIn tracking
#if In
  std::vector<std::vector<DCHitContainer> > bcInCandidates;
  BH1Filter& gFilter = BH1Filter::GetInstance();
  gFilter.Apply(*hodoAna, *DCAna, bcInCandidates);
  DCAna->TrackSearchBcIn(bcInCandidates);
  //  DCAna->TrackSearchBcIn();
  int ntIn=DCAna->GetNtracksBcIn();
  event.ntIn=ntIn;
  int ntInOk=0;
  HF1( 10, double(ntIn) );
  for( int it=0; it<ntIn; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcIn(it);
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
    
    if( chisqr < 30 ) ntInOk++;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-100;
      HF1( 13, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+1, wire-0.5 );
      HF1( 100*layerId+2, dt );
      HF1( 100*layerId+3, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+4, pos );
      HF1( 100*layerId+5, res );
      HF2( 100*layerId+6, pos, res );
      HF2( 100*layerId+7, xcal, ycal);
    }
  }
  //if( ntIn<1 ) return true;
  if( ntInOk<1 ) return true;
#endif 
  HF1( 1, 17. );

  HF1( 1, 20. );

  //if( !(ntIn==1 && ntOut ==1) ) return true;

  HF1( 1, 21. );

  //////////////K18 tracking
#if K18
  DCAna->TrackSearchK18();

  int ntK18=DCAna->GetNTracksK18();
  event.ntK18 = ntK18;
  HF1( 50, double(ntK18) );
  for( int i=0; i<ntK18; ++i ){
    K18Track *tp=DCAna->GetK18Track(i);
    if(!tp) continue;
    int nh=tp->GetNHitsTotal();
    double chisqr=tp->chisquare();
    double xt=tp->Xtgt(), yt=tp->Ytgt();
    double ut=tp->Utgt(), vt=tp->Vtgt();
    double p=tp->P();
    double cost = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*Rad2Deg;
    double phi   = atan2( ut, vt );

    HF1( 51, double(nh) );
    HF1( 52, chisqr );
    HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
    HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
    HF1( 61, p );

    event.nhK18[i] = nh;
    event.chisqrK18[i] = chisqr;
    event.p[i] = p;

    event.xtgt[i] = xt;
    event.ytgt[i] = yt;
    event.utgt[i] = ut;
    event.vtgt[i] = vt; 

    event.theta[i] = theta;
    event.phi[i] = phi;

       for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetK18HitTotal(j);
      if(!hit) continue;
      int layerId=hit->GetLayer()-100;
      HF1( 53, layerId );
      double wire=hit->GetHit()->GetWire();
      double dt=hit->GetHit()->GetDriftTime();
      double dl=hit->GetHit()->GetDriftLength();
      double pos=hit->GetCalLPos(), res=hit->GetResidual();
      DCLTrackHit *lhit=hit->GetHit();
      double xcal=lhit->GetXcal(), ycal=lhit->GetYcal();
      HF1( 100*layerId+11, double(wire)-0.5 ); 
      HF1( 100*layerId+12, dt );  HF1( 100*layerId+13, dl ); 
      HF1( 100*layerId+14, pos ); HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal );
    }

    DCLocalTrack *trIn =tp->TrackIn();
    DCLocalTrack *trOut=tp->TrackOut();
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
      double x0out=trOut->GetX(0.), y0out=trOut->GetY(0.);
      double u0out=trOut->GetU0(), v0out=trOut->GetV0();
      double chiout=trOut->GetChiSquare();
      HF1( 41, double(nhOut) ); HF1( 42, chiout );
      HF1( 44, x0out ); HF1( 45, y0out ); HF1( 46, u0out ); HF1( 47, v0out );
      HF2( 48, x0out, u0out ); HF2( 49, y0out, v0out );
      for( int jout=0; jout<nhOut; ++jout ){
	int layer=trOut->GetHit(jout)->GetLayer();
	HF1( 43, layer );
      }
    }
  }

  for( int i=0; i<ntK18; ++i ){
    K18Track *tp=DCAna->GetK18Track(i);
    if(!tp) continue;
    DCLocalTrack *trIn =tp->TrackIn();
    DCLocalTrack *trOut=tp->TrackOut();
    if( !trIn || !trOut ) continue;

    double yin=trIn->GetY(500.), vin=trIn->GetV0();
    double yout=trOut->GetY(3800.), vout=trOut->GetV0();
    HF2( 20021, yin, yout ); HF2( 20022, vin, vout );
    HF2( 20023, vin, yout ); HF2( 20024, vout, yin );
 }
#endif

  HF1( 1, 22. );

  tree->Fill();

  return true;
}

void EventK18Tracking::InitializeEvent( void )
{

  event.trigtype = -1;  
  event.nhBh2    = -1;
  event.nhBh1    = -1;
  event.ntIn     = -1;
  event.ntOut    = -1;
  event.ntK18    = -1;

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }  

  for( int it=0; it<MaxHits; it++){
    event.Bh2Seg[it] = -1;
    event.tBh2[it]  = -9999.0;
    event.deBh2[it] = -9999.0;

    event.Bh1Seg[it] = -1;
    event.tBh1[it]  = -9999.0;
    event.deBh1[it] = -9999.0;
    event.btof[it]  = -1;
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
    event.nhK18[it]     = -1; 
    event.chisqrK18[it] = -1.0;
    event.p[it]     = -9999.0;

    event.xtgt[it]  = -9999.0;
    event.ytgt[it]  = -9999.0;
    event.utgt[it]  = -9999.0;
    event.vtgt[it]  = -9999.0;
    event.theta[it] = -9999.0;
    event.phi[it]   = -9999.0;
  }
}

bool EventK18Tracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventK18Tracking;
}

const int NBin1DTBcIn  =    40;
const double MinDTBcIn =  -50.;
const double MaxDTBcIn =  250.;
const int NBin1DLBcIn  =  100;
const double MinDLBcIn = -1.0;
const double MaxDLBcIn =  3.0;

const int NBin1DTBcOut =   100;
const double MinDTBcOut = -30.;
const double MaxDTBcOut =  80.;
const int NBin1DLBcOut =   100;
const double MinDLBcOut = -1.0;
const double MaxDLBcOut =  3.0;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 30, 0., 30. );

  HB1( 10, "#Tracks BcIn", 10, 0., 10. );
  HB1( 11, "#Hits of Track BcIn", 15, 0., 15. );
  HB1( 12, "Chisqr BcIn", 500, 0., 50. ); 
  HB1( 13, "LayerId BcIn", 15, 0., 15. );
  HB1( 14, "X0 BcIn", 400, -100., 100. ); 
  HB1( 15, "Y0 BcIn", 400, -100., 100. );
  HB1( 16, "U0 BcIn", 200, -0.20, 0.20 );
  HB1( 17, "V0 BcIn", 200, -0.20, 0.20 );
  HB2( 18, "U0%X0 BcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 19, "V0%Y0 BcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 20, "X0%Y0 BcIn", 100, -100., 100., 100, -100., 100. );

  HB1( 21, "#Hits of Track BcIn [K18Track]", 15, 0., 15. );
  HB1( 22, "Chisqr BcIn [K18Track]", 500, 0., 50. ); 
  HB1( 23, "LayerId BcIn [K18Track]", 15, 0., 15. );
  HB1( 24, "X0 BcIn [K18Track]", 400, -100., 100. ); 
  HB1( 25, "Y0 BcIn [K18Track]", 400, -100., 100. );
  HB1( 26, "U0 BcIn [K18Track]", 200, -0.20, 0.20 );
  HB1( 27, "V0 BcIn [K18Track]", 200, -0.20, 0.20 );
  HB2( 28, "U0%X0 BcIn [K18Track]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 29, "V0%Y0 BcIn [K18Track]", 100, -100., 100., 100, -0.20, 0.20 );
  //HB2( 30, "X0%Y0 BcIn [K18Track]", 100, -100., 100., 100, -100., 100. );

  HB1( 30, "#Tracks BcOut", 10, 0., 10. );
  HB1( 31, "#Hits of Track BcOut", 20, 0., 20. );
  HB1( 32, "Chisqr BcOut", 500, 0., 50. ); 
  HB1( 33, "LayerId BcOut", 15, 12., 27. );
  HB1( 34, "X0 BcOut", 400, -100., 100. ); 
  HB1( 35, "Y0 BcOut", 400, -100., 100. );
  HB1( 36, "U0 BcOut",  200, -0.20, 0.20 );
  HB1( 37, "V0 BcOut",  200, -0.20, 0.20 );
  HB2( 38, "U0%X0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 39, "V0%Y0 BcOut", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 40, "X0%Y0 BcOut", 100, -100., 100., 100, -100., 100. );

  HB1( 41, "#Hits of Track BcOut [K18Track]", 20, 0., 20. );
  HB1( 42, "Chisqr BcOut [K18Track]", 500, 0., 50. ); 
  HB1( 43, "LayerId BcOut [K18Track]", 15, 12., 27. );
  HB1( 44, "X0 BcOut [K18Track]", 400, -100., 100. ); 
  HB1( 45, "Y0 BcOut [K18Track]", 400, -100., 100. );
  HB1( 46, "U0 BcOut [K18Track]",  200, -0.20, 0.20 );
  HB1( 47, "V0 BcOut [K18Track]",  200, -0.20, 0.20 );
  HB2( 48, "U0%X0 BcOut [K18Track]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 49, "V0%Y0 BcOut [K18Track]", 100, -100., 100., 100, -0.20, 0.20 );
  //HB2( 50, "X0%Y0 BcOut [K18Track]", 100, -700., 700., 100, -500., 500. );

  HB1( 50, "#Tracks K18", 10, 0., 10. );
  HB1( 51, "#Hits of K18Track", 30, 0., 30. );
  HB1( 52, "Chisqr K18Track", 500, 0., 100. );
  HB1( 53, "LayerId K18Track", 50, 0., 50. );
  HB1( 54, "Xtgt K18Track", 200, -100., 100. );
  HB1( 55, "Ytgt K18Track", 200, -100., 100. );
  HB1( 56, "Utgt K18Track", 300, -0.30, 0.30 );
  HB1( 57, "Vtgt K18Track", 300, -0.20, 0.20 );
  HB2( 58, "U%Xtgt K18Track", 100, -100., 100., 100, -0.25, 0.25 );
  HB2( 59, "V%Ytgt K18Track", 100, -100., 100., 100, -0.10, 0.10 );
  HB2( 60, "Y%Xtgt K18Track", 100, -100., 100., 100, -100., 100. );
  HB1( 61, "P K18Track", 500, 0.50, 2.0 );

  // BC1&2
  for( int i=1; i<=12; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPat Bc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 256, 0., 256. );
    title2 << "DriftTime Bc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), NBin1DTBcIn, MinDTBcIn, MaxDTBcIn );
    title3 << "DriftLength Bc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLBcIn, MinDLBcIn, MaxDLBcIn );

    std::ostringstream title4, title5, title6, title7;
    title4 << "Position Bc" << std::setw(2) << i;
    HB1( 100*i+4, title4.str().c_str(), 500, -100., 100. );
    title5 << "Residual Bc" << std::setw(2) << i;
    HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << "Resid%Pos Bc" << std::setw(2) << i;
    HB2( 100*i+6, title6.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
    title7 << "Y%Xcal Bc" << std::setw(2) << i;
    HB2( 100*i+7, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );

    title1 << " [K18Track]";
    HB1( 100*i+11, title1.str().c_str(), 256, 0., 256. );
    title2 << " [K18Track]";
    HB1( 100*i+12, title2.str().c_str(), NBin1DTBcIn, MinDTBcIn, MaxDTBcIn );
    title3 << " [K18Track]";
    HB1( 100*i+13, title3.str().c_str(), NBin1DLBcIn, MinDLBcIn, MaxDLBcIn );
    title4 << " [K18Track]";
    HB1( 100*i+14, title4.str().c_str(), 500, -100., 100. );
    title5 << " [K18Track]";
    HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << " [K18Track]";
    HB2( 100*i+16, title6.str().c_str(), 50, -100., 100., 50, -2.0, 2.0 );
    title7 << " [K18Track]";
    HB2( 100*i+17, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );
  }

  // BC34
  for( int i=13; i<=24; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "HitPat Bc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 64, 0., 64. );
    title2 << "DriftTime Bc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), NBin1DTBcOut, MinDTBcOut, MaxDTBcOut );
    title3 << "DriftLength Bc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLBcOut, MinDLBcOut, MaxDLBcOut );

    std::ostringstream title4, title5, title6, title7;
    title4 << "Position Bc" << std::setw(2) << i;
    HB1( 100*i+4, title4.str().c_str(), 1000, -100., 100. );
    title5 << "Residual Bc" << std::setw(2) << i;
    HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << "Resid%Pos Bc" << std::setw(2) << i;
    HB2( 100*i+6, title6.str().c_str(), 100, -100., 100., 50, -1.0, 1.0 );
    title7 << "Y%Xcal Bc" << std::setw(2) << i;
    HB2( 100*i+7, title7.str().c_str(), 100, -100., 100., 100, -50., 50. );

    title1 << " [K18Track]";
    HB1( 100*i+11, title1.str().c_str(), 64, 0., 64. );
    title2 << " [K18Track]";
    HB1( 100*i+12, title2.str().c_str(), NBin1DTBcOut, MinDTBcOut, MaxDTBcOut );
    title3 << " [K18Track]";
    HB1( 100*i+13, title3.str().c_str(), NBin1DLBcOut, MinDLBcOut, MaxDLBcOut );
    title4 << " [K18Track]";
    HB1( 100*i+14, title4.str().c_str(), 1000, -100., 100. );
    title5 << " [K18Track]";
    HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
    title6 << " [K18Track]";
    HB2( 100*i+16, title6.str().c_str(), 100, -100., 100., 50, -2.0, 2.0 );
    title7 << " [K18Track]";
    HB2( 100*i+17, title7.str().c_str(), 100, -100., 100., 100, -50., 50. );
  }

  HB2( 20001, "Yout%Yin", 100, -150., 150., 120, -300., 300. );
  HB2( 20002, "Vout%Vin", 100, -0.05, 0.05, 100, -0.1, 0.1 );
  HB2( 20003, "Yout%Vin", 100, -0.05, 0.05, 100, -300., 300. );
  HB2( 20004, "Yin%Vout", 100, -0.10, 0.10, 100, -150., 150. );

  HB2( 20021, "Yout%Yin [K18Track]", 100, -150., 150., 120, -300., 300. );
  HB2( 20022, "Vout%Vin [K18Track]", 100, -0.05, 0.05, 100, -0.1, 0.1 );
  HB2( 20023, "Yout%Vin [K18Track]", 100, -0.05, 0.05, 100, -300., 300. );
  HB2( 20024, "Yin%Vout [K18Track]", 100, -0.10, 0.10, 100, -150., 150. );
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");
  
 //Hodoscope
  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

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

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("p",           event.p,         "p[ntK18]/D");

  tree->Branch("xtgt",    event.xtgt,   "xtgt[ntK18]/D");
  tree->Branch("ytgt",    event.ytgt,   "ytgt[ntK18]/D");
  tree->Branch("utgt",    event.utgt,   "utgt[ntK18]/D");
  tree->Branch("vtgt",    event.vtgt,   "vtgt[ntK18]/D");

  tree->Branch("theta",   event.theta,  "theta[ntK18]/D");
  tree->Branch("phi",     event.phi,    "phi[ntK18]/D");

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

  if( K18MatrixFileName_!="" )
    K18Matrix_ = new K18TransMatrix(K18MatrixFileName_);
  if(K18Matrix_) K18Matrix_->Initialize();

  if ( bh1FilterFileName_!="" )
    BH1Filter::GetInstance().Initialize(bh1FilterFileName_);

  return true;
}
