/*
  UserBeamThrough.cc
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

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double BH2toTgt = -3.0;//Time difference between BH2 and Target [ns]

VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBeamThrough 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventBeamThrough();
  ~EventBeamThrough();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventBeamThrough::EventBeamThrough()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventBeamThrough::~EventBeamThrough()
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
  int ntInb;
  int nhInb[MaxHits]; 
  double chisqrInb[MaxHits];
  double u0Inb[MaxHits];
  double v0Inb[MaxHits];
  double x0Inb[MaxHits];
  double y0Inb[MaxHits];

  int ntOutb;
  int nhOutb[MaxHits]; 
  double chisqrOutb[MaxHits];
  double u0Outb[MaxHits];
  double v0Outb[MaxHits];
  double x0Outb[MaxHits];
  double y0Outb[MaxHits];

  int ntK18;
  int nhK18[MaxHits]; 
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double Delta[MaxHits];

  double xtgtb[MaxHits];  
  double ytgtb[MaxHits];  
  double utgtb[MaxHits];  
  double vtgtb[MaxHits]; 
  double thetab[MaxHits];
  double phib[MaxHits]; 

  //SKS
  int ntIns;
  int nhIns[MaxHits]; 
  double chisqrIns[MaxHits];
  double u0Ins[MaxHits];
  double v0Ins[MaxHits];
  double x0Ins[MaxHits];
  double y0Ins[MaxHits];

  int ntOuts;
  int nhOuts[MaxHits]; 
  double chisqrOuts[MaxHits];
  double u0Outs[MaxHits];
  double v0Outs[MaxHits];
  double x0Outs[MaxHits];
  double y0Outs[MaxHits];

  int ntSks;
  int nhSks[MaxHits]; 
  double chisqrSks[MaxHits];
  double path[MaxHits];
  double pSks[MaxHits];
  double m2[MaxHits];

  double xtgts[MaxHits];  
  double ytgts[MaxHits];  
  double utgts[MaxHits];  
  double vtgts[MaxHits]; 
  double thetas[MaxHits]; 
  double phis[MaxHits]; 

  int ndP; 
  double dP[MaxHits]; 
};
static Event event;

bool EventBeamThrough::ProcessingBegin()

{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventBeamThrough::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();
  
  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //Misc
  int trig_type=0;
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
      }//p beam
      else if( seg==2 && T>0 ){
	trig_type =(trig_type | (1 << 1));
	event.trigtype = (trig_type >> 2)+2;
      }//K baem
      else if( seg==3 && T>0 ){
	trig_type =(trig_type | (1 << 2));
	event.trigtype = (trig_type >> 3)+3;
      }//(pi, K)
      else if( seg==5 && T>0 ){
	trig_type =(trig_type | (1 << 4));
	event.trigtype = (trig_type >> 5)+5;
      }//(pi, pi)
      else if( seg==6 && T>0 ){
	trig_type =(trig_type | (1 << 5));
	event.trigtype = (trig_type >> 6)+6;
      }//(pi, p)
      else if( seg==7 && T>0 ){
	trig_type =(trig_type | (1 << 6));
	event.trigtype = (trig_type >> 7)+7;
      }
      event.trigflag[i] = T;
    }
  }

  //////////////BH2 Analysis
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  //------------------------Cut
  if( ncBh2==0 ) return true;

  //////////////BH2 time 0
//   hodoAna->DecodeBH2Hits(rawData);
//   BH2Cluster *clBH2Time0=hodoAna->GetClusterBH2(0);
//   double time0=clBH2Time0->CTime0();
//   double time0 =0.0;
//   if( clBH2Time0 ) time0=clBH2Time0->CTime0();
//   else return 0;

  //////////////TOF Number of cluster
  hodoAna->DecodeTOFHits(rawData);
  int ncTof=hodoAna->GetNClustersTOF();

  DCAna->DecodeRawHits( rawData ); 
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 

  double multi_BcIn=0.;
  double multi_BcOut=0.;
  double multi_SdcIn=0.;
  double multi_SdcOut=0.;
  
  double MaxMultiHitBcIn   = 1.0;
  double MaxMultiHitBcOut  = 100.;
  double MaxMultiHitSdcIn  = 100.;
  double MaxMultiHitSdcOut = 1.0;

  //////////////BC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetBcInHC(layer);
      int nhIn=contIn.size();
      multi_BcIn += double(nhIn);
    }
  }
  if( multi_BcIn/double(NumOfLayersBcIn) > MaxMultiHitBcIn )
   return true;
 
  //////////////BC3&4 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
    }
  }
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    return true;
 
  //////////////SDC1&2 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn += double(nhIn);
    }
  }
  if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
   return true;

  //////////////SDC3&4 number of hit in one layer not 0
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi_SdcOut += double(nhOut);
    }
  }
  if( multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut )
    return true;

  //////////////SDCIn tracking
  DCAna->TrackSearchSdcIn();

  int ntIns=DCAna->GetNtracksSdcIn();
  event.ntIns=ntIns;
  HF1( 10, double(ntIns) );
  for( int it=0; it<ntIns; ++it ){
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
    
    event.nhIns[it] = nh;
    event.chisqrIns[it] = chisqr;
    event.x0Ins[it] = x0;
    event.y0Ins[it] = y0;
    event.u0Ins[it] = u0;
    event.v0Ins[it] = v0;
   
    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer();
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
  if( ntIns<1 ) return true;

  //////////////SDCOut tracking
  DCAna->TrackSearchSdcOut();

  int ntOuts=DCAna->GetNtracksSdcOut(); 
  event.ntOuts=ntOuts;
  HF1( 30, double(ntOuts) );
  for( int it=0; it<ntOuts; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX(zTof), y0=tp->GetY(zTof);

    HF1( 31, double(nh) );
    HF1( 32, chisqr ); 
    HF1( 34, x0 ); HF1( 35, y0 );
    HF1( 36, u0 ); HF1( 37, v0 );
    HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
    HF2( 40, x0, y0 );

    event.nhOuts[it] = nh;
    event.chisqrOuts[it] = chisqr;
    event.x0Outs[it] = x0;
    event.y0Outs[it] = y0;
    event.u0Outs[it] = u0;
    event.v0Outs[it] = v0;
   
    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer();
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
  if( ntOuts<1 ) return true;

  //////////////SKS tracking
  DCAna->TrackSearchSks();

  int ntSks=DCAna->GetNTracksSks();
  event.ntSks = ntSks;
  HF1( 50, double(ntSks) );
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
    event.pSks[i] = p;

    event.xtgts[i] = xt;
    event.ytgts[i] = yt;
    event.utgts[i] = ut;
    event.vtgts[i] = vt; 

    event.thetas[i] = theta;
    event.phis[i] = phi;
    
    //std::cout<<"******************************"<<std::endl;
//     double m2;
//     for( int j=0; j<ncTof; ++j ){
//       HodoCluster *clTof=hodoAna->GetClusterTOF(j);
//       if( !clTof ) continue;
//       double time = clTof->CMeanTime()-time0-BH2toTgt;

//       //      if( time>0 ){
//       if( 1 ){
// 	HF1( 63, m2 );
// 	m2=MassSquare( p, pathL, time );
// // 		std::cout<<"Mom= "<< p <<std::endl;
// // 		std::cout<<"Path= "<< pathL <<std::endl;
// // 		std::cout<<"Time= "<< time <<std::endl;
// // 		std::cout<<"m2= "<< m2 <<std::endl;
//       }
//     }      
//     event.m2[i] = m2;
   
    for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetHit(j);
      if(!hit) continue;
      int layerId=hit->GetLayer();
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
    }
  }
  if( ntSks<1 ) return true;

//   for( int i1=0; i1<ntIns; ++i1 ){
//     DCLocalTrack *trIn=DCAna->GetTrackSdcIn(i1);
//     double yin=trIn->GetY(500.), vin=trIn->GetV0();
//     for( int i2=0; i2<ntOut; ++i2 ){
//       DCLocalTrack *trOut=DCAna->GetTrackSdcOut(i2);
//       double yout=trOut->GetY(3800.), vout=trOut->GetV0();
//       HF2( 20001, yin, yout ); HF2( 20002, vin, vout );
//       HF2( 20003, vin, yout ); HF2( 20004, vout, yin );
//     }
//   }

  for( int i=0; i<ntSks; ++i ){
    SksTrack *tp=DCAna->GetSksTrack(i);
    if(!tp) continue;
    DCLocalTrack *trIn =tp->GetLocalTrackIn();
    DCLocalTrack *trOut=tp->GetLocalTrackOut();
    if( !trIn || !trOut ) continue;

    double yin=trIn->GetY(500.), vin=trIn->GetV0();
    double yout=trOut->GetY(3800.), vout=trOut->GetV0();
    HF2( 20021, yin, yout ); HF2( 20022, vin, vout );
    HF2( 20023, vin, yout ); HF2( 20024, vout, yin );
  }

  //////////////BCOut tracking
  DCAna->TrackSearchBcOut();

  int ntOutb=DCAna->GetNtracksBcOut(); 
  event.ntOutb=ntOutb;
  HF1( 30, double(ntOutb) );
  for( int it=0; it<ntOutb; ++it ){
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

    event.nhOutb[it] = nh;
    event.chisqrOutb[it] = chisqr;
    event.x0Outb[it] = x0;
    event.y0Outb[it] = y0;
    event.u0Outb[it] = u0;
    event.v0Outb[it] = v0;
   
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
  if( ntOutb<1 ) return true;

  //////////////BCIn tracking
  std::vector<std::vector<DCHitContainer> > bcInCandidates;
  BH1Filter& gFilter = BH1Filter::GetInstance();
  gFilter.Apply(*hodoAna, *DCAna, bcInCandidates);
  DCAna->TrackSearchBcIn(bcInCandidates);
  //DCAna->TrackSearchBcIn();
  
  int ntInb=DCAna->GetNtracksBcIn();
  event.ntInb=ntInb;
  HF1( 10, double(ntInb) );
  for( int it=0; it<ntInb; ++it ){
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
    
    event.nhInb[it] = nh;
    event.chisqrInb[it] = chisqr;
    event.x0Inb[it] = x0;
    event.y0Inb[it] = y0;
    event.u0Inb[it] = u0;
    event.v0Inb[it] = v0;
   
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
  if( ntInb<1 ) return true;
  
  //////////////K18 tracking
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
    event.pK18[i] = p;

    event.xtgtb[i] = xt;
    event.ytgtb[i] = yt;
    event.utgtb[i] = ut;
    event.vtgtb[i] = vt; 

    event.thetab[i] = theta;
    event.phib[i] = phi;

       for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetK18HitTotal(j);
      if(!hit) continue;
      int layerId=hit->GetLayer();
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
  
  if( ntK18<1 ) return true;


  ///dP
  int ndP=0;
  for( int iIn=0; iIn<ntK18; ++iIn ){
    K18Track *trIn=DCAna->GetK18Track(iIn);
    if(!trIn) continue;
    double chiIn=trIn->chisquare();
    double xb=trIn->Xtgt(), yb=trIn->Ytgt();
    double ub=trIn->Utgt(), vb=trIn->Vtgt();
    double pb=trIn->P();
    double ua=trIn->TrackIn()->GetU0(), va=trIn->TrackIn()->GetV0();
    double pbp=pb-0.008-0.0200*ub;
    
    for( int iOut=0; iOut<ntSks; ++iOut ){
      SksTrack *trOut=DCAna->GetSksTrack(iOut);
      if(!trOut) continue;
      double chiOut=trOut->chisqr();
      ThreeVector Ppos=trOut->PrimaryPosition();
      ThreeVector Pmom=trOut->PrimaryMomentum();
      double ps=Pmom.mag();
      double xs=Ppos.x(), ys=Ppos.y();
      double us=Pmom.x()/Pmom.z(), vs=Pmom.y()/Pmom.z();
      double dp=ps-pb;
      double dpp=ps-pbp;


      event.dP[ndP] = dp;
      ndP++;
//       std::cout<<"***************" << std::endl;
//       std::cout<<"dP=" << dp << std::endl;
    }
  }
  event.ndP=ndP;

  tree->Fill();

  return true;
}

void EventBeamThrough::InitializeEvent( void )
{

  event.trigtype  = -1;
  event.ntInb     = -1;
  event.ntOutb    = -1;
  event.ntK18     = -1;
  event.ntIns     = -1;
  event.ntOuts    = -1;
  event.ntSks     = -1;
  event.ndP       = -1;

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }
    
  for( int it=0; it<MaxHits; it++){
    event.nhInb[it]     = -1;
    event.chisqrInb[it] = -1.0;
    event.x0Inb[it] = -9999.0;
    event.y0Inb[it] = -9999.0;
    event.u0Inb[it] = -9999.0;
    event.v0Inb[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhOutb[it]     = -1; 
    event.chisqrOutb[it] = -1.0;
    event.x0Outb[it] = -9999.0;
    event.y0Outb[it] = -9999.0;
    event.u0Outb[it] = -9999.0;
    event.v0Outb[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhIns[it]     = -1;
    event.chisqrIns[it] = -1.0;
    event.x0Ins[it] = -9999.0;
    event.y0Ins[it] = -9999.0;
    event.u0Ins[it] = -9999.0;
    event.v0Ins[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhOuts[it]     = -1; 
    event.chisqrOuts[it] = -1.0;
    event.x0Outs[it] = -9999.0;
    event.y0Outs[it] = -9999.0;
    event.u0Outs[it] = -9999.0;
    event.v0Outs[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = -1; 
    event.chisqrK18[it] = -1.0;
    event.pK18[it]     = -9999.0;
    
    event.xtgtb[it]  = -9999.0;
    event.ytgtb[it]  = -9999.0;
    event.utgtb[it]  = -9999.0;
    event.vtgtb[it]  = -9999.0;
    event.thetab[it] = -9999.0;
    event.phib[it]   = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhSks[it]     = -1; 
    event.chisqrSks[it] = -1.0;
    event.path[it]  = -9999.0;
    event.pSks[it]     = -9999.0;
    event.m2[it]    = -9999.0;

    event.xtgts[it]  = -9999.0;
    event.ytgts[it]  = -9999.0;
    event.utgts[it]  = -9999.0;
    event.vtgts[it]  = -9999.0;
    event.thetas[it] = -9999.0;
    event.phis[it]   = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.dP[it]     = -9999.0; 
  }
}

bool EventBeamThrough::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBeamThrough;
}

const int NBin1DTSdcIn  =  90;
const double MinDTSdcIn = -10.;
const double MaxDTSdcIn =  80.;
const int NBin1DLSdcIn  =  100;
const double MinDLSdcIn = -1.0;
const double MaxDLSdcIn =  3.0;

const int NBin1DTSdcOut =   400;
const double MinDTSdcOut = -100.;
const double MaxDTSdcOut =  300.;
const int NBin1DLSdcOut =   100;
const double MinDLSdcOut = -5.0;
const double MaxDLSdcOut = 15.0;

bool ConfMan:: InitializeHistograms()
{  
//   HB1( 10, "#Tracks SdcIn", 10, 0., 10. );
//   HB1( 11, "#Hits of Track SdcIn", 15, 0., 15. );
//   HB1( 12, "Chisqr SdcIn", 500, 0., 50. ); 
//   HB1( 13, "LayerId SdcIn", 15, 0., 15. );
//   HB1( 14, "X0 SdcIn", 400, -100., 100. ); 
//   HB1( 15, "Y0 SdcIn", 400, -100., 100. );
//   HB1( 16, "U0 SdcIn", 200, -0.20, 0.20 );
//   HB1( 17, "V0 SdcIn", 200, -0.20, 0.20 );
//   HB2( 18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
//   HB2( 19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
//   HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100., 100. );

//   HB1( 21, "#Hits of Track SdcIn [SksTrack]", 15, 0., 15. );
//   HB1( 22, "Chisqr SdcIn [SksTrack]", 500, 0., 50. ); 
//   HB1( 23, "LayerId SdcIn [SksTrack]", 15, 0., 15. );
//   HB1( 24, "X0 SdcIn [SksTrack]", 400, -100., 100. ); 
//   HB1( 25, "Y0 SdcIn [SksTrack]", 400, -100., 100. );
//   HB1( 26, "U0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
//   HB1( 27, "V0 SdcIn [SksTrack]", 200, -0.20, 0.20 );
//   HB2( 28, "U0%X0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
//   HB2( 29, "V0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -0.20, 0.20 );
//   //HB2( 30, "X0%Y0 SdcIn [SksTrack]", 100, -100., 100., 100, -100., 100. );

//   HB1( 30, "#Tracks SdcOut", 10, 0., 10. );
//   HB1( 31, "#Hits of Track SdcOut", 20, 0., 20. );
//   HB1( 32, "Chisqr SdcOut", 500, 0., 50. ); 
//   HB1( 33, "LayerId SdcOut", 20, 30., 50. );
//   HB1( 34, "X0 SdcOut", 1400, -1200., 1200. ); 
//   HB1( 35, "Y0 SdcOut", 1000, -500., 500. );
//   HB1( 36, "U0 SdcOut",  700, -0.35, 0.35 );
//   HB1( 37, "V0 SdcOut",  200, -0.20, 0.20 );
//   HB2( 38, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40 );
//   HB2( 39, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20 );
//   HB2( 40, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500., 500. );

//   HB1( 41, "#Hits of Track SdcOut [SksTrack]", 20, 0., 20. );
//   HB1( 42, "Chisqr SdcOut [SksTrack]", 500, 0., 50. ); 
//   HB1( 43, "LayerId SdcOut [SksTrack]", 20, 30., 50. );
//   HB1( 44, "X0 SdcOut [SksTrack]", 1400, -1200., 1200. ); 
//   HB1( 45, "Y0 SdcOut [SksTrack]", 1000, -500., 500. );
//   HB1( 46, "U0 SdcOut [SksTrack]",  700, -0.35, 0.35 );
//   HB1( 47, "V0 SdcOut [SksTrack]",  200, -0.10, 0.10 );
//   HB2( 48, "U0%X0 SdcOut [SksTrack]", 120, -600., 600., 100, -0.40, 0.40 );
//   HB2( 49, "V0%Y0 SdcOut [SksTrack]", 100, -500., 500., 100, -0.10, 0.10 );
//   //HB2( 50, "X0%Y0 SdcOut [SksTrack]", 100, -700., 700., 100, -500., 500. );

//   HB1( 50, "#Tracks SKS", 10, 0., 10. );
//   HB1( 51, "#Hits of SksTrack", 30, 0., 30. );
//   HB1( 52, "Chisqr SksTrack", 500, 0., 100. );
//   HB1( 53, "LayerId SksTrack", 50, 0., 50. );
//   HB1( 54, "Xtgt SksTrack", 200, -100., 100. );
//   HB1( 55, "Ytgt SksTrack", 200, -100., 100. );
//   HB1( 56, "Utgt SksTrack", 300, -0.30, 0.30 );
//   HB1( 57, "Vtgt SksTrack", 300, -0.20, 0.20 );
//   HB2( 58, "U%Xtgt SksTrack", 100, -100., 100., 100, -0.25, 0.25 );
//   HB2( 59, "V%Ytgt SksTrack", 100, -100., 100., 100, -0.10, 0.10 );
//   HB2( 60, "Y%Xtgt SksTrack", 100, -100., 100., 100, -100., 100. );
//   HB1( 61, "P SksTrack", 400, 0.50, 0.90 );
//   HB1( 62, "PathLength SksTrack", 600, 3000., 6000. );
//   HB1( 63, "MassSqr", 600, -0.2, 1.2 );

//   // SDC1
//   for( int i=1; i<=4; ++i ){
//     std::ostringstream title1, title2, title3;
//     title1 << "HitPat Sdc" << std::setw(2) << i;
//     HB1( 100*i+1, title1.str().c_str(), 64, 0., 64. );
//     title2 << "DriftTime Sdc" << std::setw(2) << i;
//     HB1( 100*i+2, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
//     title3 << "DriftLength Sdc" << std::setw(2) << i;
//     HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );

//     std::ostringstream title4, title5, title6, title7;
//     title4 << "Position Sdc" << std::setw(2) << i;
//     HB1( 100*i+4, title4.str().c_str(), 500, -100., 100. );
//     title5 << "Residual Sdc" << std::setw(2) << i;
//     HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << "Resid%Pos Sdc" << std::setw(2) << i;
//     HB2( 100*i+6, title6.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
//     title7 << "Y%Xcal Sdc" << std::setw(2) << i;
//     HB2( 100*i+7, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );

//     title1 << " [SksTrack]";
//     HB1( 100*i+11, title1.str().c_str(), 64, 0., 64. );
//     title2 << " [SksTrack]";
//     HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
//     title3 << " [SksTrack]";
//     HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
//     title4 << " [SksTrack]";
//     HB1( 100*i+14, title4.str().c_str(), 500, -100., 100. );
//     title5 << " [SksTrack]";
//     HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << " [SksTrack]";
//     HB2( 100*i+16, title6.str().c_str(), 50, -100., 100., 50, -2.0, 2.0 );
//     title7 << " [SksTrack]";
//     HB2( 100*i+17, title7.str().c_str(), 100, -100., 100., 50, -50., 50. );
//   }

//   // SDC2
//   for( int i=5; i<=10; ++i ){
//     std::ostringstream title1, title2, title3;
//     title1 << "HitPat Sdc" << std::setw(2) << i;
//     HB1( 100*i+1, title1.str().c_str(), 96, 0., 96. );
//     title2 << "DriftTime Sdc" << std::setw(2) << i;
//     HB1( 100*i+2, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
//     title3 << "DriftLength Sdc" << std::setw(2) << i;
//     HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );

//     std::ostringstream title4, title5, title6, title7;
//     title4 << "Position Sdc" << std::setw(2) << i;
//     HB1( 100*i+4, title4.str().c_str(), 1000, -200., 200. );
//     title5 << "Residual Sdc" << std::setw(2) << i;
//     HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << "Resid%Pos Sdc" << std::setw(2) << i;
//     HB2( 100*i+6, title6.str().c_str(), 100, -200., 200., 50, -1.0, 1.0 );
//     title7 << "Y%Xcal Sdc" << std::setw(2) << i;
//     HB2( 100*i+7, title7.str().c_str(), 100, -200., 200., 50, -100., 100. );

//     title1 << " [SksTrack]";
//     HB1( 100*i+11, title1.str().c_str(), 96, 0., 96. );
//     title2 << " [SksTrack]";
//     HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
//     title3 << " [SksTrack]";
//     HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
//     title4 << " [SksTrack]";
//     HB1( 100*i+14, title4.str().c_str(), 1000, -250., 250. );
//     title5 << " [SksTrack]";
//     HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << " [SksTrack]";
//     HB2( 100*i+16, title6.str().c_str(), 100, -200., 200., 50, -2.0, 2.0 );
//     title7 << " [SksTrack]";
//     HB2( 100*i+17, title7.str().c_str(), 100, -250., 250., 50, -100., 100. );
//   }

//   // SDC34
//   for( int i=31; i<=42; ++i ){
//     std::ostringstream title1, title2, title3;
//     title1 << "HitPat Sdc" << std::setw(2) << i;
//     HB1( 100*i+1, title1.str().c_str(), 120, 0., 120. );
//     title2 << "DriftTime Sdc" << std::setw(2) << i;
//     HB1( 100*i+2, title2.str().c_str(), NBin1DTSdcOut, MinDTSdcOut, MaxDTSdcOut );
//     title3 << "DriftLength Sdc" << std::setw(2) << i;
//     HB1( 100*i+3, title3.str().c_str(), NBin1DLSdcOut, MinDLSdcOut, MaxDLSdcOut );

//     std::ostringstream title4, title5, title6, title7;
//     title4 << "Position Sdc" << std::setw(2) << i;
//     HB1( 100*i+4, title4.str().c_str(), 1000, -600., 600. );
//     title5 << "Residual Sdc" << std::setw(2) << i;
//     HB1( 100*i+5, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << "Resid%Pos Sdc" << std::setw(2) << i;
//     HB2( 100*i+6, title6.str().c_str(), 100, -600., 600., 50, -1.0, 1.0 );
//     title7 << "Y%Xcal Sdc" << std::setw(2) << i;
//     HB2( 100*i+7, title7.str().c_str(), 100, -600., 600., 100, -600., 600. );

//     title1 << " [SksTrack]";
//     HB1( 100*i+11, title1.str().c_str(), 120, 0., 120. );
//     title2 << " [SksTrack]";
//     HB1( 100*i+12, title2.str().c_str(), NBin1DTSdcOut, MinDTSdcOut, MaxDTSdcOut );
//     title3 << " [SksTrack]";
//     HB1( 100*i+13, title3.str().c_str(), NBin1DLSdcOut, MinDLSdcOut, MaxDLSdcOut );
//     title4 << " [SksTrack]";
//     HB1( 100*i+14, title4.str().c_str(), 1000, -600., 600. );
//     title5 << " [SksTrack]";
//     HB1( 100*i+15, title5.str().c_str(), 200, -2.0, 2.0 );
//     title6 << " [SksTrack]";
//     HB2( 100*i+16, title6.str().c_str(), 100, -600., 600., 50, -2.0, 2.0 );
//     title7 << " [SksTrack]";
//     HB2( 100*i+17, title7.str().c_str(), 100, -600., 600., 100, -600., 600. );
//   }

//   HB2( 20001, "Yout%Yin", 100, -150., 150., 120, -300., 300. );
//   HB2( 20002, "Vout%Vin", 100, -0.05, 0.05, 100, -0.1, 0.1 );
//   HB2( 20003, "Yout%Vin", 100, -0.05, 0.05, 100, -300., 300. );
//   HB2( 20004, "Yin%Vout", 100, -0.10, 0.10, 100, -150., 150. );

//   HB2( 20021, "Yout%Yin [SksTrack]", 100, -150., 150., 120, -300., 300. );
//   HB2( 20022, "Vout%Vin [SksTrack]", 100, -0.05, 0.05, 100, -0.1, 0.1 );
//   HB2( 20023, "Yout%Vin [SksTrack]", 100, -0.05, 0.05, 100, -300., 300. );
//   HB2( 20024, "Yin%Vout [SksTrack]", 100, -0.10, 0.10, 100, -150., 150. );
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");

  //Beam  
  tree->Branch("ntInb",    &event.ntInb,     "ntInb/I");
  tree->Branch("nhInb",     event.nhInb,     "nhInb[ntInb]/I");
  tree->Branch("chisqrInb", event.chisqrInb, "chisqrInb[ntInb]/D");
  tree->Branch("x0Inb",     event.x0Inb,     "x0Inb[ntInb]/D");
  tree->Branch("y0Inb",     event.y0Inb,     "y0Inb[ntInb]/D");
  tree->Branch("u0Inb",     event.u0Inb,     "u0Inb[ntInb]/D");
  tree->Branch("v0Inb",     event.v0Inb,     "v0Inb[ntInb]/D");

  tree->Branch("ntOutb",   &event.ntOutb,     "ntOutb/I");
  tree->Branch("nhOutb",    event.nhOutb,     "nhOutb[ntOutb]/I");
  tree->Branch("chisqrOutb",event.chisqrOutb, "chisqrOutb[ntOutb]/D");
  tree->Branch("x0Outb",    event.x0Outb,     "x0Outb[ntOutb]/D");
  tree->Branch("y0Outb",    event.y0Outb,     "y0Outb[ntOutb]/D");
  tree->Branch("u0Outb",    event.u0Outb,     "u0Outb[ntOutb]/D");
  tree->Branch("v0Outb",    event.v0Outb,     "v0Outb[ntOutb]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");

  tree->Branch("xtgtb",    event.xtgtb,   "xtgtb[ntK18]/D");
  tree->Branch("ytgtb",    event.ytgtb,   "ytgtb[ntK18]/D");
  tree->Branch("utgtb",    event.utgtb,   "utgtb[ntK18]/D");
  tree->Branch("vtgtb",    event.vtgtb,   "vtgtb[ntK18]/D");

  tree->Branch("thetab",   event.thetab,  "thetab[ntK18]/D");
  tree->Branch("phib",     event.phib,    "phib[ntK18]/D");

  //SKS
  tree->Branch("ntIns",    &event.ntIns,     "ntIns/I");
  tree->Branch("nhIns",     event.nhIns,     "nhIns[ntIns]/I");
  tree->Branch("chisqrIns", event.chisqrIns, "chisqrIns[ntIns]/D");
  tree->Branch("x0Ins",     event.x0Ins,     "x0Ins[ntIns]/D");
  tree->Branch("y0Ins",     event.y0Ins,     "y0Ins[ntIns]/D");
  tree->Branch("u0Ins",     event.u0Ins,     "u0Ins[ntIns]/D");
  tree->Branch("v0Ins",     event.v0Ins,     "v0Ins[ntIns]/D");

  tree->Branch("ntOuts",   &event.ntOuts,     "ntOuts/I");
  tree->Branch("nhOuts",    event.nhOuts,     "nhOuts[ntOuts]/I");
  tree->Branch("chisqrOuts",event.chisqrOuts, "chisqrOuts[ntOuts]/D");
  tree->Branch("x0Outs",    event.x0Outs,     "x0Outs[ntOuts]/D");
  tree->Branch("y0Outs",    event.y0Outs,     "y0Outs[ntOuts]/D");
  tree->Branch("u0Outs",    event.u0Outs,     "u0Outs[ntOuts]/D");
  tree->Branch("v0Outs",    event.v0Outs,     "v0Outs[ntOuts]/D");


  tree->Branch("ntSks",      &event.ntSks,     "ntSks/I");
  tree->Branch("nhSks",       event.nhSks,     "nhSks[ntSks]/I");
  tree->Branch("chisqrSks",   event.chisqrSks, "chisqrSks[ntSks]/D");
  tree->Branch("path",        event.path,      "path[ntSks]/D");
  tree->Branch("pSks",        event.pSks,      "pSks[ntSks]/D");
  tree->Branch("m2",          event.m2,        "m2[ntSks]/D");

  tree->Branch("xtgts",    event.xtgts,   "xtgts[ntSks]/D");
  tree->Branch("ytgts",    event.ytgts,   "ytgts[ntSks]/D");
  tree->Branch("utgts",    event.utgts,   "utgts[ntSks]/D");
  tree->Branch("vtgts",    event.vtgts,   "vtgts[ntSks]/D");

  tree->Branch("thetas",   event.thetas,  "thetas[ntSks]/D");
  tree->Branch("phis",     event.phis,    "phis[ntSks]/D");

  tree->Branch("ndP",    &event.ndP,    "ndP/I");
  tree->Branch("dP",      event.dP,    "dP[ndP]/D");

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
