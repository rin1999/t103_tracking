/*
  UserDCMonitor.cc
  2009/11  K.Shirotori 
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
#include "EvDisp.hh"
#include "SksLib.hh"

const double TdcLow  =  700.;
const double TdcHigh = 1000.;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);


 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventEvDisp 
  : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
  HodoAnalyzer *hodoAna;
public:
  EventEvDisp();
  ~EventEvDisp();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
};

EventEvDisp::EventEvDisp()
  : VEvent(),
    rawData(0), DCAna(new DCAnalyzer()), hodoAna(new HodoAnalyzer)
{
}

EventEvDisp::~EventEvDisp()
{
  if (DCAna)   delete DCAna;
  if (hodoAna) delete hodoAna;
  if (rawData) delete rawData;
}

bool EventEvDisp::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventEvDisp::ProcessingNormal()
{
  bool FlagEvDisp=false;
  FlagEvDisp=ConfMan::GetConfManager()->GetEvDispFlag();



  rawData = new RawData;
  rawData->DecodeHits();

  //BH2
  int bh2_nhits=0;
  {
    const HodoRHitContainer &cont=rawData->GetBH2RawHC();
    int nh=cont.size();
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId();
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if (Tu>0 || Td>0) {
        if (FlagEvDisp) {
          const EvDisp & evDisp = EvDisp::GetInstance();          
          evDisp.DrawHitBH2(seg, Tu, Td);
        }
      }
    }
  }

  //TOF
  int tof_nhits=0;
  {
    int TofId=DCGeomMan::GetInstance().GetTofId(); 
    const HodoRHitContainer &cont=rawData->GetTOFRawHC();
    int nh=cont.size();
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId();
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if (Tu>0 || Td>0) {
        if (FlagEvDisp) {
          const EvDisp & evDisp = EvDisp::GetInstance();          
          evDisp.DrawHitHodoscope(TofId, seg, Tu, Td);
        }
      }
    }
  }

  //AC
  int ac1_nhits=0;
  int ac2_nhits=0;
  {  
    for( int layer=0; layer<NumOfLayersAc; ++layer ){
      const HodoRHitContainer &cont=rawData->GetACRawHC(layer);
      int nh=cont.size();
      int nh1=0, nh2=0;
      for( int i=0; i<nh; ++i ){
	HodoRawHit *hit=cont[i];
	int seg=hit->SegmentId();
	int A=hit->GetAdc1();
	int T=hit->GetTdc1();
      }
    }
  }

  //LC
  int lc_nhits=0;
  {
    int LcId=DCGeomMan::GetInstance().GetLcId();
    const HodoRHitContainer &cont=rawData->GetLCRawHC();
    int nh=cont.size();
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId();
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if (Tu>0 || Td>0) {
	if (FlagEvDisp) {
	  const EvDisp & evDisp = EvDisp::GetInstance();
	  evDisp.DrawHitHodoscope(LcId, seg, Tu, Td);
	}
      }
    }
  }


  int multi[NumOfLayersSdcOut+1];
  for (int i=0; i<12; i++)
    multi[i]=0;

  DCAna->DecodeRawHits( rawData );
  //**************************************************************************
  //******************RawData

  //SDC3&SDC4
  /*
  {
    for( int layer=0; layer<NumOfLayersSdcOut; ++layer ){
      if( layer<NumOfLayersSdc ){
	const DCRHitContainer &contOut =rawData->GetSdcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 100*(layer+1), nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+1)+1, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    HF1( 100*(layer+1)+2, tdc );
	    HF1( 10000*(layer+1)+wire, tdc );
	  }
	}
      }
      if( layer>NumOfLayersSdc-1 ){
	const DCRHitContainer &contOut =rawData->GetSdcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 100*(layer+1)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+1)+1, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    HF1( 100*(layer+1)+2, tdc );
	    HF1( 10000*(layer+1)+wire, tdc );
	  }
	}
      }
    }
  }
  */

  //SDC1&SDC2

  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      for( int i=0; i<nhIn; ++i ){
	DCHit *hit=contIn[i];
	double wire=hit->GetWire();
	int multi = hit->GetTdcSize();
	bool goodFlag=false;
	for (int j=0; j<multi; j++) {
	  if (hit->rangecheck(j)) {
	    goodFlag=true;
	    break;
	  }
	}
	if (FlagEvDisp) {
	  const EvDisp & evDisp = EvDisp::GetInstance();
	  if (goodFlag)
	    evDisp.DrawHitWire(layer, int(wire));
	}

      }
    }
  }

  /*
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCRHitContainer &contIn =rawData->GetSdcInRawHC(layer);
      int nhIn=contIn.size();
      for( int i=0; i<nhIn; ++i ){
	DCRawHit *hit=contIn[i];
	int wire=hit->WireId();
	std::cout << layer << ", " << wire << std::endl;
	if (FlagEvDisp) {
	  const EvDisp & evDisp = EvDisp::GetInstance();
	  evDisp.DrawHitWire(layer, wire);
	}
      }
    }
  }
  */
  /*
  {
    for( int layer=0; layer<NumOfLayersSdcIn; ++layer ){
      const DCHitContainer &contIn =DCAna->GetSdcInHC(layer);
      int nhIn=contIn.size();
      for( int i=0; i<nhIn; ++i ){
	DCHit *hit=contIn[i];
	int wire=hit->GetWire();
	std::cout << layer << ", " << wire << std::endl;
	if (FlagEvDisp) {
	  const EvDisp & evDisp = EvDisp::GetInstance();
	  evDisp.DrawHitWire(layer+1, wire);
	}
      }
    }
  }
  */

  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      multi[layer] = nhOut;
      HF1( 100*(layer+1), nhOut );
      for( int i=0; i<nhOut; ++i ){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	
	if (FlagEvDisp) {
	  const EvDisp & evDisp = EvDisp::GetInstance();
	  evDisp.DrawHitWire(layer+30, int(wire));
	}

	HF1( 100*layer+1, wire-0.5 );
	int nhtdc = hit->GetTdcSize();
	for( int k=0; k<nhtdc; k++ ){
	  int tdc = hit->GetTdcVal(k);
	  HF1( 100*layer+2, tdc );
	  HF1( 10000*layer+int(wire), tdc );
	}
	int nhdt = hit->GetDriftTimeSize();
	for( int k=0; k<nhdt; k++ ){
	  double dt = hit->GetDriftTime(k);
	  HF1( 100*layer+3, dt );
	  HF1( 10000*layer+1000+int(wire), dt );
	}
	int nhdl = hit->GetDriftLengthSize();
	for( int k=0; k<nhdl; k++ ){
	  double dl = hit->GetDriftLength(k);
	  HF1( 100*layer+4, dl );
	  HF1( 10000*layer+2000+int(wire), dl );
	}
      }
    }
  }


  DCAna->TrackSearchSdcIn();
  int ntSdcIn=DCAna->GetNtracksSdcIn();
  for( int it=0; it<ntSdcIn; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;
    //evDisp
    if (FlagEvDisp) {
      const EvDisp & evDisp = EvDisp::GetInstance();
      evDisp.DrawSdcInLocalTrack(tp);
    }
  }


#if 0
  // SDC3
  //std::cout << "=============TrackSearchSdc3==========" << std::endl;
  DCAna->TrackSearchSdc3();
  
  int ntSdc3=DCAna->GetNtracksSdc3();
  for( int it=0; it<ntSdc3; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdc3(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-PlOffsSdcOut;
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

      //      std::cout<<"layerId = "<<layerId<<std::endl;

       if (fabs(dl-fabs(xlcal-wp))<2.0) {      
	 HFProf( 100*layerId+20, dt, fabs(xlcal-wp));
	 HF2( 100*layerId+22, dt, fabs(xlcal-wp));
	 HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	 HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
       }


       /*
       if (nt != 1)
	 continue;

       const DCHitContainer & cont = DCAna->GetBDHC(layerId-1);
       int nhits=cont.size();
       if (nhits>=2) {
	 bool flagConsistent=false;
	 for (int j=0; j<nhits; j++) {
	   DCHit *tmphit=cont[j];
	   int lid = tmphit->GetLayer();
	   int wid = tmphit->GetWire();
	   if (lid == layerId && wid == wire)
	     flagConsistent=true;
	 }

	 if (flagConsistent) {
	   for (int j=0; j<nhits; j++) {
	     DCHit *tmphit=cont[j];
	     int wid = tmphit->GetWire();
	     double dt_notTrack = tmphit->GetDriftTime();
	     if ( wid != wire) {
	       HF2( 100*layerId+52, dt, dt_notTrack);	       
	       HF1( 100*layerId+53, dt - dt_notTrack);	       

	       if (nhits == 2) {
		 if (fabs(wid-wire)==1.0) {
		   HF2( 100*layerId+54, dt, dt_notTrack);	       
		   HF1( 100*layerId+55, dt - dt_notTrack);	       
		 } else {
		   HF2( 100*layerId+56, dt, dt_notTrack);	       
		   HF1( 100*layerId+57, dt - dt_notTrack);	       
		 }
	       }
	     }
	   }
	 } else {
	   std::cout << "Hit wire is not consistent, layerId = " 
		     << layerId << ", wire = " << wire << std::endl;
	 }
	 
       }
       */
    }
  }

  // SDC4
  //std::cout << "==========TrackSearchSdc4============" << std::endl;
  DCAna->TrackSearchSdc4();
  
  int ntSdc4=DCAna->GetNtracksSdc4();
  for( int it=0; it<ntSdc4; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdc4(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-PlOffsSdcOut;
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

      //      std::cout<<"layerId = "<<layerId<<std::endl;

       if (fabs(dl-fabs(xlcal-wp))<2.0) {      
	 HFProf( 100*layerId+20, dt, fabs(xlcal-wp));
	 HF2( 100*layerId+22, dt, fabs(xlcal-wp));
	 HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	 HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
       }

       /*
       if (nt != 1)
	 continue;

       const DCHitContainer & cont = DCAna->GetBDHC(layerId-1);
       int nhits=cont.size();
       if (nhits>=2) {
	 bool flagConsistent=false;
	 for (int j=0; j<nhits; j++) {
	   DCHit *tmphit=cont[j];
	   int lid = tmphit->GetLayer();
	   int wid = tmphit->GetWire();
	   if (lid == layerId && wid == wire)
	     flagConsistent=true;
	 }

	 if (flagConsistent) {
	   for (int j=0; j<nhits; j++) {
	     DCHit *tmphit=cont[j];
	     int wid = tmphit->GetWire();
	     double dt_notTrack = tmphit->GetDriftTime();
	     if ( wid != wire) {
	       HF2( 100*layerId+52, dt, dt_notTrack);	       
	       HF1( 100*layerId+53, dt - dt_notTrack);	       

	       if (nhits == 2) {
		 if (fabs(wid-wire)==1.0) {
		   HF2( 100*layerId+54, dt, dt_notTrack);	       
		   HF1( 100*layerId+55, dt - dt_notTrack);	       
		 } else {
		   HF2( 100*layerId+56, dt, dt_notTrack);	       
		   HF1( 100*layerId+57, dt - dt_notTrack);	       
		 }
	       }
	     }
	   }
	 } else {
	   std::cout << "Hit wire is not consistent, layerId = " 
		     << layerId << ", wire = " << wire << std::endl;
	 }
	 
       }
       */
    }
  }
#endif

  for (int i=1; i<=NumOfLayersSdcOut; i++) {
    if (i==6)
      continue;
    //    if (multi[i] != 1) {
      if (FlagEvDisp) {
	const EvDisp & evDisp = EvDisp::GetInstance();
	evDisp.get_command();
	evDisp.EndOfEvent();
      }
//            return true;
//     }
  }


#if 1
  // Sdc Out
  std::cout << "==========TrackSearch SdcOut============" << std::endl;
  DCAna->TrackSearchSdcOut();
  
  int nt=DCAna->GetNtracksSdcOut();
  for( int it=0; it<nt; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    //evDisp
    if (FlagEvDisp) {
      const EvDisp & evDisp = EvDisp::GetInstance();
      evDisp.DrawSdcOutLocalTrack(tp);
    }

    double cost = 1./sqrt(1.+u0*u0+v0*v0);
    double theta = acos(cost)*Rad2Deg;

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-PlOffsSdcOut;
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

      //      std::cout<<"layerId = "<<layerId<<std::endl;

       if (fabs(dl-fabs(xlcal-wp))<2.0) {      
	 HFProf( 100*layerId+20, dt, fabs(xlcal-wp));
	 HF2( 100*layerId+22, dt, fabs(xlcal-wp));
	 HFProf( 100000*layerId+3000+int(wire), xlcal-wp,dt);
	 HF2( 100000*layerId+4000+int(wire), xlcal-wp,dt);
       }

       /*
       if (nt != 1)
	 continue;

       const DCHitContainer & cont = DCAna->GetBDHC(layerId-1);
       int nhits=cont.size();
       if (nhits>=2) {
	 bool flagConsistent=false;
	 for (int j=0; j<nhits; j++) {
	   DCHit *tmphit=cont[j];
	   int lid = tmphit->GetLayer();
	   int wid = tmphit->GetWire();
	   if (lid == layerId && wid == wire)
	     flagConsistent=true;
	 }

	 if (flagConsistent) {
	   for (int j=0; j<nhits; j++) {
	     DCHit *tmphit=cont[j];
	     int wid = tmphit->GetWire();
	     double dt_notTrack = tmphit->GetDriftTime();
	     if ( wid != wire) {
	       HF2( 100*layerId+52, dt, dt_notTrack);	       
	       HF1( 100*layerId+53, dt - dt_notTrack);	       

	       if (nhits == 2) {
		 if (fabs(wid-wire)==1.0) {
		   HF2( 100*layerId+54, dt, dt_notTrack);	       
		   HF1( 100*layerId+55, dt - dt_notTrack);	       
		 } else {
		   HF2( 100*layerId+56, dt, dt_notTrack);	       
		   HF1( 100*layerId+57, dt - dt_notTrack);	       
		 }
	       }
	     }
	   }
	 } else {
	   std::cout << "Hit wire is not consistent, layerId = " 
		     << layerId << ", wire = " << wire << std::endl;
	 }
	 
       }
       */
    }
  }
#endif

  DCAna->TrackSearchSks();

  int ntSks=DCAna->GetNTracksSks();
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
    
  }

  if (FlagEvDisp) {
    const EvDisp & evDisp = EvDisp::GetInstance();
    evDisp.get_command();
    evDisp.EndOfEvent();
  }

  return true;
}

bool EventEvDisp::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventEvDisp;
}


const int NbinSdcOutTdc   = 2000;
const double MinSdcOutTdc =    0.;
const double MaxSdcOutTdc = 2000.;

const int NbinSdcOutDT   = 200;
const double MinSdcOutDT =  -50.;
const double MaxSdcOutDT =  350.;

const int NbinSdcOutDL   = 200;
const double MinSdcOutDL =  -5.;
const double MaxSdcOutDL =  15.;

bool ConfMan:: InitializeHistograms()
{  

  //***********************Chamber
  // SDC3
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC3#" << std::setw(2) << i;
    title2 << "Hitpat SDC3#" << std::setw(2) << i;
    title3 << "Tdc SDC3#" << std::setw(2) << i;
    title4 << "Drift Time SDC3#" << std::setw(2) << i;
    title5 << "Drift Length SDC3#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 100*i+2, title3.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 100*i+3, title4.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*i+4, title5.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    int NumOfWires;
    if (i==2 || i==5)
      NumOfWires = MaxWireSDC3X;
    else 
      NumOfWires = MaxWireSDC3V;
    for (int wire=0; wire<NumOfWires; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
      title12 << "Drift Time SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title13 << "Drift Length SDC3#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );

    }

  }

  // SDC4
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC4#" << std::setw(2) << i;
    title2 << "Hitpat SDC4#" << std::setw(2) << i;
    title3 << "Tdc SDC4#" << std::setw(2) << i;
    title4 << "Drift Time SDC4#" << std::setw(2) << i;
    title5 << "Drift Length SDC4#" << std::setw(2) << i;
    HB1( 100*(i+6)+0, title1.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 100*(i+6)+3, title4.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
    HB1( 100*(i+6)+4, title5.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );
    int NumOfWires;
    if (i==2 || i==5)
      NumOfWires = MaxWireSDC4X;
    else 
      NumOfWires = MaxWireSDC4V;
    for (int wire=0; wire<NumOfWires; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+wire, title11.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
      title12 << "Drift Time SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+1000+wire, title12.str().c_str(), NbinSdcOutDT, MinSdcOutDT, MaxSdcOutDT );
      title13 << "Drift Length SDC4#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*(i+6)+2000+wire, title13.str().c_str(), NbinSdcOutDL, MinSdcOutDL, MaxSdcOutDL );

    }
  }

  // Tracking Histgrams
  for( int i=1; i<=12; ++i ){
    std::ostringstream title0, title1, title2, title3, title4;
    std::ostringstream title5, title6, title7, title8, title9;
    std::ostringstream title10, title11, title12, title13, title14;
    
    title0 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+11, title0.str().c_str(), 120, 0., 120. );
    title1 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+12, title1.str().c_str(), 500, -100, 400 );
    title2 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    HB1( 100*i+13, title2.str().c_str(), 100, -5, 15);
    title3 << "Position SdcOut" << std::setw(2) << i;
    HB1( 100*i+14, title3.str().c_str(), 100, -1000., 1000. ); 
    title4 << "Residual SdcOut" << std::setw(2) << i;
    HB1( 100*i+15, title4.str().c_str(), 200, -5.0, 5.0 );
    title5 << "Resid%Pos SdcOut" << std::setw(2) << i;
    HB2( 100*i+16, title5.str().c_str(), 400, -1000., 1000., 100, -1.0, 1.0 );
    title6 << "Y%Xcal SdcOut" << std::setw(2) << i;
    HB2( 100*i+17, title6.str().c_str(), 100, -1000., 1000., 100, -1000., 1000. );
    title7 << "Res%dl SdcOut" << std::setw(2) << i;
    HB2( 100*i+18, title7.str().c_str(), 100, -12., 12., 100, -3.0, 3.0 );
    title8 << "Hit Pos%Drift Time SdcOut" << std::setw(2) << i;
    HB2( 100*i+19, title8.str().c_str(), 100, -100., 400., 100, -20, 20);

    title9 << "Drift Length%Drift Time SdcOut" << std::setw(2) << i;
    HBProf( 100*i+20, title9.str().c_str(), 100, -50, 300, 0, 12);
    HB2( 100*i+22, title9.str().c_str(), 100, -50, 300, 100,0, 12);

    title4 << " w/o self plane"; 
    HB1( 100*i+21, title4.str().c_str(), 200, -5.0, 5.0 );

    title10 << "Residual SdcOut (0<theta<15)" << std::setw(2) << i;
    HB1( 100*i+71, title10.str().c_str(), 200, -5.0, 5.0 );

    title11 << "Residual SdcOut (15<theta<30)" << std::setw(2) << i;
    HB1( 100*i+72, title11.str().c_str(), 200, -5.0, 5.0 );

    title12 << "Residual SdcOut (30<theta<45)" << std::setw(2) << i;
    HB1( 100*i+73, title12.str().c_str(), 200, -5.0, 5.0 );

    title13 << "Residual SdcOut (45<theta)" << std::setw(2) << i;
    HB1( 100*i+74, title13.str().c_str(), 200, -5.0, 5.0 );

    for (int j=0; j<120; j++) {
      std::ostringstream title;
      title << "XT of Layer " << std::setw(2) << i << " Wire #"<< std::setw(4) << j;
      HBProf( 100000*i+3000+j, title.str().c_str(), 101, -12.12, 12.12, -30,300);
      HB2( 100000*i+4000+j, title.str().c_str(), 100, -12, 12, 100, -30,300);
    }

    /*
    title9 << " w/o self plane"; 
    HB2( 100*i+23, title9.str().c_str(), 100, -5., 5., 150, -50, 100);
    HBProf( 100*i+24, title9.str().c_str(), 100, -50, 150, 0, 5);
    HB2( 100*i+25, title9.str().c_str(), 100, -50, 150, 100,0, 3.5);
    title3 << " w/o self plane"; 
    HB2( 100*i+26, title3.str().c_str(), 50, -100., 100., 50, -1.0, 1.0 );
    */

  }
  
  return true;

}

bool ConfMan::InitializeParameterFiles( void )
{
//   if( CMapFileName_!="" )
//     CMapManager_ = new CMapMan(CMapFileName_);
//   if(CMapManager_) CMapManager_->Initialize();

//   if( HodoParamFileName_!="" )
//     HodoParamManager_ = new HodoParamMan(HodoParamFileName_);
//   if(HodoParamManager_) HodoParamManager_->Initialize();

//   if( HodoPHCFileName_!="" )
//     HodoPHCManager_ = new HodoPHCMan(HodoPHCFileName_);
//   if( HodoPHCManager_ ) HodoPHCManager_->Initialize();

//   if( ScalerDefinitionFileName_!="" )
//     ScalerAnalyzer_ = new ScalerAna(ScalerDefinitionFileName_);
//   else
//     ScalerAnalyzer_ = new ScalerAna();
//   if(ScalerAnalyzer_) ScalerAnalyzer_->Initialize(); 

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
