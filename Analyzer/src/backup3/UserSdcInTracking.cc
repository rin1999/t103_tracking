/*
  UserSdcInTracking.cc
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
const int PIK_TRIG   = 2;
const int BT_TRIG    = 3;
const int PIPI_TRIG  = 11;
const int PIP_TRIG   = 12;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxMultiHitSdcIn = 100.0;
const double MaxMultiHitBcOut = 10.0;

const double TdcLow  =    0.;
const double TdcHigh = 1100.;

#define HodoCut  1
#define BcOutCut 0

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventSdcInTracking 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;

public:
  EventSdcInTracking();
  ~EventSdcInTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventSdcInTracking::EventSdcInTracking()
  : VEvent(),
    rawData(0), 
    DCAna(new DCAnalyzer()),
    hodoAna(new HodoAnalyzer())
{
}

EventSdcInTracking::~EventSdcInTracking()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

bool EventSdcInTracking::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

#ifndef MaxHits 
#define MaxHits 20
#endif

//For Tree
struct Event{
  int trigpat[NumOfMisc];
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

  int ntrack;
  double chisqr[MaxHits];  
  double x0[MaxHits];  
  double y0[MaxHits];  
  double u0[MaxHits];  
  double v0[MaxHits];  
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

bool EventSdcInTracking::ProcessingNormal()
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
      if( seg==1 && T>0 ){ trig_type = seg; }
      //(pi, K)
      //if( seg==2 && T>0 ){ trig_type = seg;}
      //pi beam x TOF
      //if( seg==3 && T>0 ){ trig_type = seg; }
      //(pi, pi)
      //if( seg==11 && T>0 ){ trig_type = seg; }
      event.trigflag[seg] = T;      

      if( T>0 ) event.trigpat[i] = seg;
    }
  }

  HF1( 1, 0. );

  //------------------------Cut
  if (trig_type != BEAM_TRIG)  return true;
  //if (trig_type != PIK_TRIG)  return true;
  //if (trig_type != BT_TRIG)  return true;
  //if(trig_type != PIPI_TRIG)  return true;

  HF1( 1, 1. );

#if HodoCut
 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int ncBh2=hodoAna->GetNClustersBH2();
  event.nhBh2=ncBh2;
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
#endif

  HF1( 1, 9. );
  DCAna->DecodeRawHits( rawData );

#if BcOutCut
  //BC3&BC4
  double multi_BcOut=0.;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetBcOutHC(layer);
      int nhOut=contOut.size();
      multi_BcOut += double(nhOut);
      HF1( 100*layer, nhOut );
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	
	HF1( 100*layer+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	}
	HF1( 100*layer+2, tdc1st );
	HF1( 10000*layer+int(wire), tdc1st );

	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+int(wire), dt );
	}
	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+4, dl );
	}
      }
    }    
  }
  if( multi_BcOut/double(NumOfLayersBcOut) > MaxMultiHitBcOut )
    return true;

  {
    // Bc Out
    //  std::cout << "==========TrackSearch BcOut============" << std::endl;
    int ntOk = 0;
    DCAna->TrackSearchBcOut();
    int nt=DCAna->GetNtracksBcOut(); 
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();
      
      double cost = 1./sqrt(1.+u0*u0+v0*v0);
      double theta = acos(cost)*Rad2Deg;
      double xtgt=tp->GetX(1800.), ytgt=tp->GetY(1800.);
      double utgt=u0, vtgt=v0;
      double xbac=tp->GetX(885.), ybac=tp->GetY(885.);
      double ubac=u0, vbac=v0;
      
      double xsdcIn=tp->GetX(2100.), ysdcIn=tp->GetY(2100.);
      double usdcIn=u0, vsdcIn=v0;
      
      if( -80.<xsdcIn && xsdcIn<80. &&
	  -30. <ysdcIn && ysdcIn<30. && chisqr<100 )
	ntOk++;
    }
    if( ntOk<1 ) return true;
  }
#endif

  //////////////SDC1&2 number of hit layer
  HF1( 1, 10. );

  //SDC1&SDC2
  double multi_SdcIn=0.;
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      multi_SdcIn += double(nhIn);
      HF1( 100*layer, nhIn );
      int tdchits=0;
      for( int i=0; i<nhIn; ++i ){
	DCHit *hit=contIn[i];
	double wire=hit->GetWire();
	HF1( 100*layer+1, wire-0.5 );

	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ){
	    tdc1st=tdc;
	    tdchits++;
	  }
	}
	HF1( 100*layer+2, tdc1st );
	HF1( 10000*layer+int(wire), tdc1st );
	HF1( 100*layer+9, tdchits );
	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+int(wire), dt );
	}

	int nhdl = hit->GetDriftTimeSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+4, dl );
	}
      }
    }    
  }
 if( multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn )
   return true;

 HF1( 1, 11. );

#if 1
  // Sdc In
  //std::cout << "==========TrackSearch SdcIn============" << std::endl;
  DCAna->TrackSearchSdcIn();
  int nt=DCAna->GetNtracksSdcIn();
  event.ntrack=nt;
  HF1( 10, double(nt) );
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1( 11, double(nh) );
    HF1( 12, chisqr ); 
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    double xtgt=tp->GetX(0), ytgt=tp->GetY(0);
    double utgt=u0, vtgt=v0;
    HF1( 21, xtgt ); HF1( 22, ytgt );
    HF1( 23, utgt ); HF1( 24, vtgt );
    HF2( 25, xtgt, utgt ); HF2( 26, ytgt, vtgt );
    HF2( 27, xtgt, ytgt );

    double xbac=tp->GetX(-815.0), ybac=tp->GetY(-815.0);
    double ubac=u0, vbac=v0;
    HF1( 31, xbac ); HF1( 32, ybac );
    HF1( 33, ubac ); HF1( 34, vbac );
    HF2( 35, xbac, ubac ); HF2( 36, ybac, vbac );
    HF2( 37, xbac, ybac );

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer(); 
      HF1( 13, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 100*layerId+11, wire-0.5 );
      HF1( 100*layerId+12, dt );
      HF1( 100*layerId+13, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 100*layerId+14, pos );
      HF1( 100*layerId+15, res );
      HF2( 100*layerId+16, pos, res );
      HF2( 100*layerId+17, xcal, ycal);
      //      HF1( 100000*layerId+50000+wire, res);      
      double wp=hit->GetWirePosition();
      double sign=1.;
      if( pos-wp<0. ) sign=-1;
      HF2( 100*layerId+18, sign*dl, res );
      double xlcal=hit->GetLocalCalPos();
      HF2( 100*layerId+19, dt, xlcal-wp);

      if (theta>=0 && theta<15) 
	HF1( 100*layerId+71, res );
      else if (theta>=15 && theta<30) 
	HF1( 100*layerId+72, res );
      else if (theta>=30 && theta<45) 
	HF1( 100*layerId+73, res );
      else if (theta>=45) 
	HF1( 100*layerId+74, res );
      
      if (fabs(dl-fabs(xlcal-wp))<2.0) {      
	HFProf( 100*layerId+20, dt, fabs(xlcal-wp));
	HF2( 100*layerId+22, dt, fabs(xlcal-wp));
	HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
      }
      
//       if (nt != 1)
// 	continue;
      
//       const DCHitContainer & cont = DCAna->GetBDHC(layerId-1);
//       int nhits=cont.size();
//       if (nhits>=2) {
// 	bool flagConsistent=false;
// 	for (int j=0; j<nhits; j++) {
// 	  DCHit *tmphit=cont[j];
// 	  int lid = tmphit->GetLayer();
// 	  int wid = tmphit->GetWire();
// 	  if (lid == layerId && wid == wire)
// 	    flagConsistent=true;
// 	}
	
// 	if (flagConsistent) {
// 	  for (int j=0; j<nhits; j++) {
// 	    DCHit *tmphit=cont[j];
// 	    int wid = tmphit->GetWire();
// 	    double dt_notTrack = tmphit->GetDriftTime();
// 	    if ( wid != wire) {
// 	      HF2( 100*layerId+52, dt, dt_notTrack);	       
// 	      HF1( 100*layerId+53, dt - dt_notTrack);	       
	      
// 	      if (nhits == 2) {
// 		if (fabs(wid-wire)==1.0) {
// 		  HF2( 100*layerId+54, dt, dt_notTrack);	       
// 		  HF1( 100*layerId+55, dt - dt_notTrack);	       
// 		} else {
// 		  HF2( 100*layerId+56, dt, dt_notTrack);	       
// 		  HF1( 100*layerId+57, dt - dt_notTrack);	       
// 		}
// 	      }
// 	    }
// 	  }
// 	} else {
// 	  std::cout << "Hit wire is not consistent, layerId = " 
// 		    << layerId << ", wire = " << wire << std::endl;
// 	}
//       }
    }
  }
#endif

  HF1( 1, 12. );

  tree->Fill();

  return true;
}

void EventSdcInTracking::InitializeEvent( void )
{
  event.ntrack   = -1;
  event.nhBh2  = -1;
  event.nhBh1  = -1;

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
    event.chisqr[it] = -1.0;
    event.x0[it] = -9999.0;
    event.y0[it] = -9999.0;
    event.u0[it] = -9999.0;
    event.v0[it] = -9999.0;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }
}

bool EventSdcInTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventSdcInTracking;
}


const int NbinSdcInTdc   = 2000;
const double MinSdcInTdc =    0.;
const double MaxSdcInTdc = 2000.;

const int NbinSdcInDT   =   100;
const double MinSdcInDT =  -30.;
const double MaxSdcInDT =   80.;

const int NbinSdcInDL   = 100;
const double MinSdcInDL =  -1.;
const double MaxSdcInDL =   3.;

bool ConfMan:: InitializeHistograms()
{  
  HB1( 1, "Status", 20, 0., 20. );

  //***********************Chamber
  // SDC1
  for( int i=1; i<=NumOfLayersSdc-2; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    std::ostringstream title9;
    title1 << "#Hits SDC1#" << std::setw(2) << i;
    title2 << "Hitpat SDC1#" << std::setw(2) << i;
    title3 << "Tdc SDC1#" << std::setw(2) << i;
    title4 << "Drift Time SDC1#" << std::setw(2) << i;
    title5 << "Drift Length SDC1#" << std::setw(2) << i;
    title9 << "#Hits w/ TDC cut SDC1#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinSdcInDT, MinSdcInDT, MaxSdcInDT );
    HB1( 100*i+4, title5.str().c_str(), NbinSdcInDL, MinSdcInDL, MaxSdcInDL );
    HB1( 100*i+9, title9.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );

    for (int wire=1; wire<=MaxWireSDC1; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
      title12 << "Drift Time SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinSdcInDT, MinSdcInDT, MaxSdcInDT );
      title13 << "Drift Length SDC1#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinSdcInDL, MinSdcInDL, MaxSdcInDL );
    }
  }

  // SDC2
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    std::ostringstream title9;
    title1 << "#Hits SDC2#" << std::setw(2) << i;
    title2 << "Hitpat SDC2#" << std::setw(2) << i;
    title3 << "Tdc SDC2#" << std::setw(2) << i;
    title4 << "Drift Time SDC2#" << std::setw(2) << i;
    title5 << "Drift Length SDC2#" << std::setw(2) << i;
    title9 << "#Hits w/ TDC cut SDC2#" << std::setw(2) << i;

    HB1( 100*(i+4)+0, title1.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 100*(i+4)+1, title2.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 100*(i+4)+2, title3.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 100*(i+4)+3, title4.str().c_str(), NbinSdcInDT, MinSdcInDT, MaxSdcInDT );
    HB1( 100*(i+4)+4, title5.str().c_str(), NbinSdcInDL, MinSdcInDL, MaxSdcInDL );
    HB1( 100*(i+4)+9, title9.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );

    for (int wire=1; wire<=MaxWireSDC2; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+wire, title11.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
      title12 << "Drift Time SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+1000+wire, title12.str().c_str(), NbinSdcInDT, MinSdcInDT, MaxSdcInDT );
      title13 << "Drift Length SDC2#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+4)+2000+wire, title13.str().c_str(), NbinSdcInDL, MinSdcInDL, MaxSdcInDL );

    }
  }

  // Tracking Histgrams
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
  HB2( 20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100, 100 );

  HB1( 21, "Xtgt SdcIn", 400, -100., 100. ); 
  HB1( 22, "Ytgt SdcIn", 400, -100., 100. );
  HB1( 23, "Utgt SdcIn", 200, -0.20, 0.20 );
  HB1( 24, "Vtgt SdcIn", 200, -0.20, 0.20 );
  HB2( 25, "Utgt%Xtgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 26, "Vtgt%Ytgt SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 27, "Xtgt%Ytgt SdcIn", 100, -100., 100., 100, -100, 100 );

  HB1( 31, "Xbac SdcIn", 400, -100., 100. ); 
  HB1( 32, "Ybac SdcIn", 400, -100., 100. );
  HB1( 33, "Ubac SdcIn", 200, -0.20, 0.20 );
  HB1( 34, "Vbac SdcIn", 200, -0.20, 0.20 );
  HB2( 35, "Ubac%Xbac SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 36, "Vbac%Ybac SdcIn", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 37, "Xbac%Ybac SdcIn", 100, -100., 100., 100, -100, 100 );

  for( int i=1; i<=NumOfLayersSdcIn; ++i ){
    std::ostringstream title0, title1, title2, title3, title4;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 100, 0., 100. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), 200, -10, 150 );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), 50, -0.5, 3);
    title3 << "Position SdcIn" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -250., 250. ); 
    title4 << "Residual SdcIn" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 400, -2.0, 2.0 );
    title5 << "Resid%Pos SdcIn" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 250, -250., 250., 100, -1.0, 1.0 );
    title6 << "Y%Xcal SdcIn" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -250., 250., 100, -250., 250. );
    title7 << "Res%dl SdcIn" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 100, -3., 3., 100, -1.0, 1.0 );
    title8 << "Hit Pos%Drift Time SdcIn" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), 100, -5., 100., 100, -3, 3);

    title9 << "Drift Length%Drift Time SdcIn" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 100, -5, 50, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 100, -5, 50, 50, 0, 12);

    title4 << " w/o self plane"; 
    HB1( 100*i+21, title4.str().c_str(), 200, -5.0, 5.0 );

    title10 << "Residual SdcIn (0<theta<15)" << std::setw(2) << i;
    HB1( 100*i+71, title10.str().c_str(), 200, -5.0, 5.0 );

    title11 << "Residual SdcIn (15<theta<30)" << std::setw(2) << i;
    HB1( 100*i+72, title11.str().c_str(), 200, -5.0, 5.0 );

    title12 << "Residual SdcIn (30<theta<45)" << std::setw(2) << i;
    HB1( 100*i+73, title12.str().c_str(), 200, -5.0, 5.0 );

    title13 << "Residual SdcIn (45<theta)" << std::setw(2) << i;
    HB1( 100*i+74, title13.str().c_str(), 200, -5.0, 5.0 );

    for (int j=1; j<=64; j++) {
      std::ostringstream title;
      title << "XT of Layer " << std::setw(2) << i << " Wire #"<< std::setw(4) << j;
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -12.12, 12.12, -30,300);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -12, 12, 100, -30,300);
    }
    
    //     title9 << " w/o self plane"; 
    //     HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    //     HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    //     HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    //     title3 << " w/o self plane"; 
    //     HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
  }
  
  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,   "trigflag[20]/I");

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

  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chosqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

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


  return true;
}
