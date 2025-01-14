/*
  UserMonitor.cc
  2009/11 K.Shirotori 
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
#include "HodoRawHit.hh"
#include "DCRawHit.hh"

static const int GCHid =   9000;
static const int BH1Hid = 10000;
static const int BH2Hid = 20000;
static const int BACHid = 30000;
static const int TGTHid = 40000;
static const int TOFHid = 50000;
static const int ACHid  = 60000;
static const int LCHid  = 70000;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventMonitor();
  ~EventMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
};

EventMonitor::EventMonitor()
  : VEvent(),
    rawData(0)
{
}

EventMonitor::~EventMonitor()
{
  if (rawData) delete rawData;
}

bool EventMonitor::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

const double AdcMin  = 0.;
const double AdcMax  = 4096.;
const double TdcMin  = 0.;
const double TdcMax  = 4096.;

bool EventMonitor::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();

  //**************************************************************************
  //******************RawData

  //GC
  {
    const HodoRHitContainer &cont=rawData->GetGCRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
       HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int A=hit->GetAdc1();
      int T=hit->GetTdc1();

      //Up
      HF1( GCHid+100*seg+1, double(A) );
      if( T>0 ){
	HF1( GCHid+100*seg+3, double(T) );
	HF1( GCHid+100*seg+5, double(A) );
      }
      else{
	HF1( GCHid+100*seg+7, double(A) );
      }
    }
  }

  //BH1
  {
    const HodoRHitContainer &cont=rawData->GetBH1RawHC();
    int nh=cont.size();
    HF1( BH1Hid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="BH1";
      HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( BH1Hid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( BH1Hid+100*seg+3, double(Tu) );
	HF1( BH1Hid+100*seg+5, double(Au) );
      }
      else{
	HF1( BH1Hid+100*seg+7, double(Au) );
      }
      //Down
      HF1( BH1Hid+100*seg+2, double(Ad) );
      if( Tu>0 ){
	HF1( BH1Hid+100*seg+4, double(Td) );
	HF1( BH1Hid+100*seg+6, double(Ad) );
      }
      else{
	HF1( BH1Hid+100*seg+7, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( BH1Hid+3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( BH1Hid+5, seg-0.5 );
      }
    }
    HF1( BH1Hid+2, double(nh1) ); HF1( BH1Hid+4, double(nh2) );
  }

  //BH2
  {
    const HodoRHitContainer &cont=rawData->GetBH2RawHC();
    int nh=cont.size();
    HF1( BH2Hid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="BH2";
      HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( BH2Hid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( BH2Hid+100*seg+3, double(Tu) );
	HF1( BH2Hid+100*seg+5, double(Au) );
      }
      else{
	HF1( BH2Hid+100*seg+7, double(Au) );
      }
      //Down
      HF1( BH2Hid+100*seg+2, double(Ad) );
      if( Tu>0 ){
	HF1( BH2Hid+100*seg+4, double(Td) );
	HF1( BH2Hid+100*seg+6, double(Ad) );
      }
      else{
	HF1( BH2Hid+100*seg+7, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( BH2Hid+3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( BH2Hid+5, seg-0.5 );
      }
    }
    HF1( BH2Hid+2, double(nh1) ); HF1( BH2Hid+4, double(nh2) );
  }

  //BAC
  {
    const HodoRHitContainer &cont=rawData->GetBACRawHC();
    int nh=cont.size();
    HF1( BACHid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="BAC";
	HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( BACHid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( BACHid+100*seg+3, double(Tu) );
	HF1( BACHid+100*seg+5, double(Au) );
      }
      else{
	HF1( BACHid+100*seg+7, double(Au) );
      }
      //       //Down
      //       HF1( BACHid+100*seg+2, double(Ad) );
      //       if( Tu>0 ){
      // 	HF1( BACHid+100*seg+4, double(Td) );
      // 	HF1( BACHid+100*seg+6, double(Ad) );
      //       }
      //       else{
      // 	HF1( BACHid+100*seg+7, double(Ad) );
      //       }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( BACHid+3, seg-0.5 );
      }
      //       if( Tu>0 && Td>0 ){
      // 	++nh2; HF1( BACHid+5, seg-0.5 );
      //       }
    }
    HF1( BACHid+2, double(nh1) ); //HF1( BACHid+4, double(nh2) );
  }

  //TGT
  {
    const HodoRHitContainer &cont=rawData->GetTGTRawHC();
    int nh=cont.size();
    HF1( TGTHid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="TGT";
	HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( TGTHid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( TGTHid+100*seg+3, double(Tu) );
	HF1( TGTHid+100*seg+5, double(Au) );
      }
      else{
	HF1( TGTHid+100*seg+7, double(Au) );
      }
      //       //Down
      //       HF1( TGTHid+100*seg+2, double(Ad) );
      //       if( Tu>0 ){
      // 	HF1( TGTHid+100*seg+4, double(Td) );
      // 	HF1( TGTHid+100*seg+6, double(Ad) );
      //       }
      //       else{
      // 	HF1( TGTHid+100*seg+7, double(Ad) );
      //       }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( TGTHid+3, seg-0.5 );
      }
      //       if( Tu>0 && Td>0 ){
      // 	++nh2; HF1( TGTHid+5, seg-0.5 );
      //       }
    }
    HF1( TGTHid+2, double(nh1) ); //HF1( TGTHid+4, double(nh2) );
  }

  //TOF
  {
    const HodoRHitContainer &cont=rawData->GetTOFRawHC();
    int nh=cont.size();
    HF1( TOFHid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="TOF";
      HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( TOFHid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( TOFHid+100*seg+3, double(Tu) );
	HF1( TOFHid+100*seg+5, double(Au) );
      }
      else{
	HF1( TOFHid+100*seg+7, double(Au) );
      }
      //Down
      HF1( TOFHid+100*seg+2, double(Ad) );
      if( Tu>0 ){
	HF1( TOFHid+100*seg+4, double(Td) );
	HF1( TOFHid+100*seg+6, double(Ad) );
      }
      else{
	HF1( TOFHid+100*seg+7, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( TOFHid+3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( TOFHid+5, seg-0.5 );
      }
    }
    HF1( TOFHid+2, double(nh1) ); HF1( TOFHid+4, double(nh2) );
  }

  //AC
  {  
    for( int layer=0; layer<NumOfLayersAc; ++layer ){
      const HodoRHitContainer &cont=rawData->GetACRawHC(layer);
      int nh=cont.size();
      if( layer==0 ) HF1( ACHid, double(nh) );
      if( layer==1 ) HF1( ACHid, double(nh)+20. );
      int nh1=0, nh2=0;
      for( int i=0; i<nh; ++i ){
	HodoRawHit *hit=cont[i];
	int seg=hit->SegmentId()+1;
	std::string name="AC";
	HF1( name.c_str(), seg-0.5 );
	int A=hit->GetAdc1();
	int T=hit->GetTdc1();

	if( layer==0 ){
	  //AC1
	  HF1( ACHid+100*seg+1, double(A) );
	  if( T>0 ){
	    HF1( ACHid+100*seg+3, double(T) );
	    HF1( ACHid+100*seg+5, double(A) );
	  }
	  else{
	    HF1( ACHid+100*seg+7, double(A) );
	  }
	  if( T>0 ){
	    ++nh1; HF1( ACHid+3, seg-0.5 );
	  }
	  if( T>0 ){
	    ++nh2; HF1( ACHid+5, seg-0.5 );
	  }
	}
	if( layer==1 ){
	  //AC2
	  HF1( ACHid+100*seg+2, double(A) );
	  if( T>0 ){
	    HF1( ACHid+100*seg+4, double(T) );
	    HF1( ACHid+100*seg+6, double(A) );
	  }
	  else{
	    HF1( ACHid+100*seg+7, double(A) );
	  }
	  //HitPat
	  if( T>0 ){
	    ++nh1; HF1( ACHid+3, seg-0.5+20.0 );
	  }
	  if( T>0 ){
	    ++nh2; HF1( ACHid+5, seg-0.5+20.0 );
	  }
	}
      }
      if( layer==0 ){
	HF1( ACHid+2, double(nh1) ); 
	HF1( ACHid+4, double(nh2) );
      }
      if( layer==1 ){
	HF1( ACHid+2, double(nh1)+20.0 ); 
	HF1( ACHid+4, double(nh2)+20.0 );
      }
    }
  }

  //LC
  {
    const HodoRHitContainer &cont=rawData->GetLCRawHC();
    int nh=cont.size();
    HF1( LCHid, double(nh) );
    int nh1=0, nh2=0;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="LC";
      HF1( name.c_str(), seg-0.5 );
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Up
      HF1( LCHid+100*seg+1, double(Au) );
      if( Tu>0 ){
	HF1( LCHid+100*seg+3, double(Tu) );
	HF1( LCHid+100*seg+5, double(Au) );
      }
      else{
	HF1( LCHid+100*seg+7, double(Au) );
      }
      //Down
      HF1( LCHid+100*seg+2, double(Ad) );
      if( Tu>0 ){
	HF1( LCHid+100*seg+4, double(Td) );
	HF1( LCHid+100*seg+6, double(Ad) );
      }
      else{
	HF1( LCHid+100*seg+7, double(Ad) );
      }
      //HitPat
      if( Tu>0 || Td>0 ){
	++nh1; HF1( LCHid+3, seg-0.5 );
      }
      if( Tu>0 && Td>0 ){
	++nh2; HF1( LCHid+5, seg-0.5 );
      }
    }
    HF1( LCHid+2, double(nh1) ); HF1( LCHid+4, double(nh2) );
  }

  //SDC3&SDC4
  {
    for( int layer=0; layer<NumOfLayersSdcOut; ++layer ){
      if( layer<NumOfLayersSdc ){
	const DCRHitContainer &contIn =rawData->GetSdcOutRawHC(layer);
	int nhIn=contIn.size();
	HF1( 3000+100*(layer+1+10)+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+1+10)+2, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( tdc>tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*(layer+1+10)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*(layer+1+10)+1, wire-0.5 );
	  HF1( 3000+100*(layer+1+10)+4, tdc1st );
	}
      }
      if( layer>NumOfLayersSdc-1 ){
	const DCRHitContainer &contIn =rawData->GetSdcOutRawHC(layer);
	int nhIn=contIn.size();
	HF1( 3000+100*(layer+1+16)+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+1+10)+2, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( tdc > tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*(layer+1+16)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*(layer+1+16)+1, wire-0.5 );
	  HF1( 3000+100*(layer+1+16)+4, tdc1st );
	}
      }
    }
  }

  return true;
}

bool EventMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor;
}

const int NbinAdc = 1024;
const double MinAdc  =    0.;
const double MaxAdc  = 4096.;

const int NbinTdc = 1024;
const double MinTdc  =    0.;
const double MaxTdc  = 4096.;

const int NbinBcInTdc   =  1024;
const double MinBcInTdc =    0.;
const double MaxBcInTdc = 4096.;

const int NbinBcOutTdc   =  1024;
const double MinBcOutTdc =    0.;
const double MaxBcOutTdc = 4096.;

const int NbinSdcInTdc   =  1024;
const double MinSdcInTdc =    0.;
const double MaxSdcInTdc = 4096.;

const int NbinSdcOutTdc   =  500;
const double MinSdcOutTdc =    0.;
const double MaxSdcOutTdc = 2000.;

bool ConfMan:: InitializeHistograms()
{  
  //***********************HodoScope
  // GC
  HB1( GCHid, "#Hits GC",        NumOfSegGC+1, 0., double(NumOfSegGC+1) );
  HB1( GCHid+1, "Hitpat GC",       NumOfSegGC,   0., double(NumOfSegGC)   );
  HB1( GCHid+2, "#Hits GC(Tor)",   NumOfSegGC+1, 0., double(NumOfSegGC+1) );
  HB1( GCHid+3, "Hitpat GC(Tor)",  NumOfSegGC,   0., double(NumOfSegGC)   );
  HB1( GCHid+4, "#Hits GC(Tand)",  NumOfSegGC+1, 0., double(NumOfSegGC+1) );
  HB1( GCHid+5, "Hitpat GC(Tand)", NumOfSegGC,   0., double(NumOfSegGC)   );

  for( int i=1; i<=NumOfSegGC; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "GC-" << i << " UpAdc";
    HB1( GCHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "GC-" << i << " DownAdc";
    HB1( GCHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "GC-" << i << " UpTdc";
    HB1( GCHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "GC-" << i << " DownTdc";
    HB1( GCHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "GC-" << i << " UpAdc(w Tdc)";
    HB1( GCHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "GC-" << i << " DownAdc(w Tdc)";
    HB1( GCHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    std::stringstream title7, title8;
    title7 << "GC-" << i << " UpAdc(w/o Tdc)";
    HB1( GCHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "GC-" << i << " DownAdc(w/o Tdc)";
    HB1( GCHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // BH1
  HB1( BH1Hid, "#Hits BH1",        NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+1, "Hitpat BH1",       NumOfSegBH1,   0., double(NumOfSegBH1)   );
  HB1( BH1Hid+2, "#Hits BH1(Tor)",   NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+3, "Hitpat BH1(Tor)",  NumOfSegBH1,   0., double(NumOfSegBH1)   );
  HB1( BH1Hid+4, "#Hits BH1(Tand)",  NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+5, "Hitpat BH1(Tand)", NumOfSegBH1,   0., double(NumOfSegBH1)   );

  for( int i=1; i<=NumOfSegBH1; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BH1-" << i << " UpAdc";
    HB1( BH1Hid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "BH1-" << i << " DownAdc";
    HB1( BH1Hid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "BH1-" << i << " UpTdc";
    HB1( BH1Hid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "BH1-" << i << " DownTdc";
    HB1( BH1Hid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "BH1-" << i << " UpAdc(w Tdc)";
    HB1( BH1Hid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "BH1-" << i << " DownAdc(w Tdc)";
    HB1( BH1Hid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "BH1-" << i << " UpAdc(w/o Tdc)";
    HB1( BH1Hid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "BH1-" << i << " DownAdc(w/o Tdc)";
    HB1( BH1Hid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // BH2
  HB1( BH2Hid, "#Hits BH2",        NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+1, "Hitpat BH2",       NumOfSegBH2,   0., double(NumOfSegBH2)   );
  HB1( BH2Hid+2, "#Hits BH2(Tor)",   NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+3, "Hitpat BH2(Tor)",  NumOfSegBH2,   0., double(NumOfSegBH2)   );
  HB1( BH2Hid+4, "#Hits BH2(Tand)",  NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+5, "Hitpat BH2(Tand)", NumOfSegBH2,   0., double(NumOfSegBH2)   );

  for( int i=1; i<=NumOfSegBH2; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BH2-" << i << " UpAdc";
    HB1( BH2Hid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "BH2-" << i << " DownAdc";
    HB1( BH2Hid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "BH2-" << i << " UpTdc";
    HB1( BH2Hid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "BH2-" << i << " DownTdc";
    HB1( BH2Hid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "BH2-" << i << " UpAdc(w Tdc)";
    HB1( BH2Hid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "BH2-" << i << " DownAdc(w Tdc)";
    HB1( BH2Hid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "BH2-" << i << " UpAdc(w/o Tdc)";
    HB1( BH2Hid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "BH2-" << i << " DownAdc(w/o Tdc)";
    HB1( BH2Hid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // BAC
  HB1( BACHid, "#Hits BAC",        NumOfSegBAC+1, 0., double(NumOfSegBAC+1) );
  HB1( BACHid+1, "Hitpat BAC",       NumOfSegBAC,   0., double(NumOfSegBAC)   );
  HB1( BACHid+2, "#Hits BAC(Tor)",   NumOfSegBAC+1, 0., double(NumOfSegBAC+1) );
  HB1( BACHid+3, "Hitpat BAC(Tor)",  NumOfSegBAC,   0., double(NumOfSegBAC)   );
  HB1( BACHid+4, "#Hits BAC(Tand)",  NumOfSegBAC+1, 0., double(NumOfSegBAC+1) );
  HB1( BACHid+5, "Hitpat BAC(Tand)", NumOfSegBAC,   0., double(NumOfSegBAC)   );

  for( int i=1; i<=NumOfSegBAC; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BAC-" << i << " UpAdc";
    HB1( BACHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "BAC-" << i << " DownAdc";
    HB1( BACHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "BAC-" << i << " UpTdc";
    HB1( BACHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "BAC-" << i << " DownTdc";
    HB1( BACHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "BAC-" << i << " UpAdc(w Tdc)";
    HB1( BACHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "BAC-" << i << " DownAdc(w Tdc)";
    HB1( BACHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "BAC-" << i << " UpAdc(w/o Tdc)";
    HB1( BACHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "BAC-" << i << " DownAdc(w/o Tdc)";
    HB1( BACHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // TGT
  HB1( TGTHid, "#Hits TGT",        NumOfSegTGT+1, 0., double(NumOfSegTGT+1) );
  HB1( TGTHid+1, "Hitpat TGT",       NumOfSegTGT,   0., double(NumOfSegTGT)   );
  HB1( TGTHid+2, "#Hits TGT(Tor)",   NumOfSegTGT+1, 0., double(NumOfSegTGT+1) );
  HB1( TGTHid+3, "Hitpat TGT(Tor)",  NumOfSegTGT,   0., double(NumOfSegTGT)   );
  HB1( TGTHid+4, "#Hits TGT(Tand)",  NumOfSegTGT+1, 0., double(NumOfSegTGT+1) );
  HB1( TGTHid+5, "Hitpat TGT(Tand)", NumOfSegTGT,   0., double(NumOfSegTGT)   );

  for( int i=1; i<=NumOfSegTGT; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "TGT-" << i << " UpAdc";
    HB1( TGTHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "TGT-" << i << " DownAdc";
    HB1( TGTHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "TGT-" << i << " UpTdc";
    HB1( TGTHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "TGT-" << i << " DownTdc";
    HB1( TGTHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "TGT-" << i << " UpAdc(w Tdc)";
    HB1( TGTHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "TGT-" << i << " DownAdc(w Tdc)";
    HB1( TGTHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "TGT-" << i << " UpAdc(w/o Tdc)";
    HB1( TGTHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "TGT-" << i << " DownAdc(w/o Tdc)";
    HB1( TGTHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // TOF
  HB1( TOFHid, "#Hits TOF",        NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+1, "Hitpat TOF",       NumOfSegTOF,   0., double(NumOfSegTOF)   );
  HB1( TOFHid+2, "#Hits TOF(Tor)",   NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+3, "Hitpat TOF(Tor)",  NumOfSegTOF,   0., double(NumOfSegTOF)   );
  HB1( TOFHid+4, "#Hits TOF(Tand)",  NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+5, "Hitpat TOF(Tand)", NumOfSegTOF,   0., double(NumOfSegTOF)   );

  for( int i=1; i<=NumOfSegTOF; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "TOF-" << i << " UpAdc";
    HB1( TOFHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "TOF-" << i << " DownAdc";
    HB1( TOFHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "TOF-" << i << " UpTdc";
    HB1( TOFHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "TOF-" << i << " DownTdc";
    HB1( TOFHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "TOF-" << i << " UpAdc(w Tdc)";
    HB1( TOFHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "TOF-" << i << " DownAdc(w Tdc)";
    HB1( TOFHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "TOF-" << i << " UpAdc(w/o Tdc)";
    HB1( TOFHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "TOF-" << i << " DownAdc(w/o Tdc)";
    HB1( TOFHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // AC
  HB1( ACHid, "#Hits AC",        2*NumOfSegAC+1, 0., double(2*NumOfSegAC+1) );
  HB1( ACHid+1, "Hitpat AC",       2*NumOfSegAC,   0., double(2*NumOfSegAC)   );
  HB1( ACHid+2, "#Hits AC(Tor)",   2*NumOfSegAC+1, 0., double(2*NumOfSegAC+1) );
  HB1( ACHid+3, "Hitpat AC(Tor)",  2*NumOfSegAC,   0., double(2*NumOfSegAC)   );
  HB1( ACHid+4, "#Hits AC(Tand)",  2*NumOfSegAC+1, 0., double(2*NumOfSegAC+1) );
  HB1( ACHid+5, "Hitpat AC(Tand)", 2*NumOfSegAC,   0., double(2*NumOfSegAC)   );

  for( int i=1; i<=2*NumOfSegAC; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "AC-" << i << " UpAdc";
    HB1( ACHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "AC-" << i << " DownAdc";
    HB1( ACHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "AC-" << i << " UpTdc";
    HB1( ACHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "AC-" << i << " DownTdc";
    HB1( ACHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "AC-" << i << " UpAdc(w Tdc)";
    HB1( ACHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "AC-" << i << " DownAdc(w Tdc)";
    HB1( ACHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "AC-" << i << " UpAdc(w/o Tdc)";
    HB1( ACHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "AC-" << i << " DownAdc(w/o Tdc)";
    HB1( ACHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  // LC
  HB1( LCHid, "#Hits LC",        NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+1, "Hitpat LC",       NumOfSegLC,   0., double(NumOfSegLC)   );
  HB1( LCHid+2, "#Hits LC(Tor)",   NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+3, "Hitpat LC(Tor)",  NumOfSegLC,   0., double(NumOfSegLC)   );
  HB1( LCHid+4, "#Hits LC(Tand)",  NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+5, "Hitpat LC(Tand)", NumOfSegLC,   0., double(NumOfSegLC)   );

  for( int i=1; i<=NumOfSegLC; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "LC-" << i << " UpAdc";
    HB1( LCHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "LC-" << i << " DownAdc";
    HB1( LCHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "LC-" << i << " UpTdc";
    HB1(  LCHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "LC-" << i << " DownTdc";
    HB1( LCHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "LC-" << i << " UpAdc(w Tdc)";
    HB1( LCHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "LC-" << i << " DownAdc(w Tdc)";
    HB1( LCHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "LC-" << i << " UpAdc(w/o Tdc)";
    HB1( LCHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "LC-" << i << " DownAdc(w/o Tdc)";
    HB1( LCHid+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
  }

  //***********************Chamber
  // BC1
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits BC1#" << std::setw(2) << i;
    title2 << "Hitpat BC1#" << std::setw(2) << i;
    title3 << "Tdc BC1#" << std::setw(2) << i;
    HB1( 100*(i+0), title1.str().c_str(), MaxWireBC1+1, 0., double(MaxWireBC1+1) );
    HB1( 100*(i+0)+1, title2.str().c_str(), MaxWireBC1+1, 0., double(MaxWireBC1+1) );
    HB1( 100*(i+0)+2, title3.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
  }

  // BC2
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits BC2#" << std::setw(2) << i;
    title2 << "Hitpat BC2#" << std::setw(2) << i;
    title3 << "Tdc BC2#" << std::setw(2) << i;
    HB1( 100*(i+6), title1.str().c_str(), MaxWireBC2+1, 0., double(MaxWireBC2+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireBC2+1, 0., double(MaxWireBC2+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
  }

  // BC3
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits BC3#" << std::setw(2) << i;
    title2 << "Hitpat BC3#" << std::setw(2) << i;
    title3 << "Tdc BC3#" << std::setw(2) << i;
    HB1( 100*(i+12), title1.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*(i+12)+1, title2.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*(i+12)+2, title3.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
  }

  // BC4
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits BC4#" << std::setw(2) << i;
    title2 << "Hitpat BC4#" << std::setw(2) << i;
    title3 << "Tdc BC4#" << std::setw(2) << i;
    HB1( 100*(i+18), title1.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+18)+1, title2.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+18)+2, title3.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
  }
  
  // SDC1
  for( int i=1; i<=NumOfLayersSdc-2; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits SDC1#" << std::setw(2) << i;
    title2 << "Hitpat SDC1#" << std::setw(2) << i;
    title3 << "Tdc SDC1#" << std::setw(2) << i;
    HB1( 3000+100*(i+0), title1.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 3000+100*(i+0)+1, title2.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 3000+100*(i+0)+2, title3.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
  }

  // SDC2
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3;
    title1 << "#Hits SDC2#" << std::setw(2) << i;
    title2 << "Hitpat SDC2#" << std::setw(2) << i;
    title3 << "Tdc SDC2#" << std::setw(2) << i;
    HB1( 3000+100*(i+4), title1.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 3000+100*(i+4)+1, title2.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 3000+100*(i+4)+2, title3.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
  }

  // SDC3
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC3#" << std::setw(2) << i;
    title2 << "Hitpat SDC3#" << std::setw(2) << i;
    title3 << "TDC First Hitpat SDC3#" << std::setw(2) << i;
    title4 << "Tdc SDC3#" << std::setw(2) << i;
    title5 << "Tdc First SDC3#" << std::setw(2) << i;
    HB1( 3000+100*(i+10)+0, title1.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 3000+100*(i+10)+1, title2.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );
    HB1( 3000+100*(i+10)+2, title3.str().c_str(), MaxWireSDC3+1, 0., double(MaxWireSDC3+1) );

    HB1( 3000+100*(i+10)+3, title4.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 3000+100*(i+10)+4, title5.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
  }

  // SDC4
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC4#" << std::setw(2) << i;
    title2 << "Hitpat SDC4#" << std::setw(2) << i;
    title3 << "TDC First Hitpat SDC4#" << std::setw(2) << i;
    title4 << "Tdc SDC4#" << std::setw(2) << i;
    title5 << "Tdc First SDC4#" << std::setw(2) << i;
    HB1( 3000+100*(i+16)+0, title1.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 3000+100*(i+16)+1, title2.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );
    HB1( 3000+100*(i+16)+2, title3.str().c_str(), MaxWireSDC4+1, 0., double(MaxWireSDC4+1) );

    HB1( 3000+100*(i+16)+3, title4.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
    HB1( 3000+100*(i+16)+4, title5.str().c_str(), NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc );
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

  return true;
}
