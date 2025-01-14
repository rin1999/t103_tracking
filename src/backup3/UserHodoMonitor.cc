/*
  UserHodoMonitor.cc
  2009/11  K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoRawHit.hh"

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

class EventHodoMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventHodoMonitor();
  ~EventHodoMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  void InitializeEvent();
  bool InitializeHistograms();
};

EventHodoMonitor::EventHodoMonitor()
  : VEvent(),
    rawData(0)
{
}

EventHodoMonitor::~EventHodoMonitor()
{
  if (rawData) delete rawData;
}

bool EventHodoMonitor::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

#ifndef MaxHits 
#define MaxHits 64
#endif

//For Tree
struct Event{
  int trigpat[MaxHits];  
  int trigflag[NumOfMisc];

  int gcnhits;
  int gchitpat[MaxHits];
  double gca[NumOfSegGC];
  double gct[NumOfSegGC];

  int bh1nhits;
  int bh1hitpat[MaxHits];
  double bh1ua[NumOfSegBH1];
  double bh1ut[NumOfSegBH1];
  double bh1da[NumOfSegBH1];
  double bh1dt[NumOfSegBH1];

  int bh2nhits;
  int bh2hitpat[MaxHits];
  double bh2ua[NumOfSegBH2];
  double bh2ut[NumOfSegBH2];
  double bh2da[NumOfSegBH2];
  double bh2dt[NumOfSegBH2];

  int bacnhits;
  int bachitpat[MaxHits];
  double bacua[NumOfSegBAC];
  double bacut[NumOfSegBAC];
  double bacda[NumOfSegBAC];
  double bacdt[NumOfSegBAC];

  int tgtnhits;
  int tgthitpat[MaxHits];
  double tgta[NumOfSegTGT];
  double tgtt[NumOfSegTGT];

  int tofnhits;
  int tofhitpat[MaxHits];
  double tofua[NumOfSegTOF];
  double tofut[NumOfSegTOF];
  double tofda[NumOfSegTOF];
  double tofdt[NumOfSegTOF];

  int ac1nhits;
  int ac2nhits;
  int ac1hitpat[MaxHits];
  int ac2hitpat[MaxHits];
  double ac1a[NumOfSegAC];
  double ac1t[NumOfSegAC];
  double ac2a[NumOfSegAC];
  double ac2t[NumOfSegAC];

  int lcnhits;
  int lchitpat[MaxHits];
  double lcua[NumOfSegLC];
  double lcut[NumOfSegLC];
  double lcda[NumOfSegLC];
  double lcdt[NumOfSegLC];
};
static Event event;

const double AdcMin  = 0.;
const double AdcMax  = 4096.;
const double TdcMin  = 0.;
const double TdcMax  = 4096.;

bool EventHodoMonitor::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();

  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //**************************************************************************
  //******************RawData

  //Misc
  {
    const HodoRHitContainer &cont=rawData->GetMiscRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int T=hit->GetTdc1();
      
      if( T>0 ) event.trigpat[i] = seg;  
      event.trigflag[i] = T;
    }
  }

  //GC
  int gc_nhits=0;
  {
    const HodoRHitContainer &cont=rawData->GetGCRawHC();
    int nh=cont.size();
    HF1( GCHid, double(nh) );
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      int A=hit->GetAdc1();
      int T=hit->GetTdc1();
  
      //Tree
      if( T>0 ){
	event.gchitpat[i] = seg;
	gc_nhits++; 
      }
      event.gca[i] = A;
      event.gct[i] = T;

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
  event.gcnhits = gc_nhits;

  //BH1
  int bh1_nhits=0;
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

      //Tree
      if( Tu>0 || Td>0 ){
	event.bh1hitpat[i] = seg;
	bh1_nhits++; 
      }
      event.bh1ua[i] = Au;
      event.bh1ut[i] = Tu;
      event.bh1da[i] = Ad;
      event.bh1dt[i] = Td;

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
      if( Td>0 ){
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
  event.bh1nhits = bh1_nhits;

  //BH2
  int bh2_nhits=0;
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

      //Tree
      if( Tu>0 || Td>0 ){
	event.bh2hitpat[i]= seg;
	bh2_nhits++; 
      }
      event.bh2ua[i] = Au;
      event.bh2ut[i] = Tu;
      event.bh2da[i] = Ad;
      event.bh2dt[i] = Td;

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
      if( Td>0 ){
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
  event.bh2nhits =bh2_nhits;

  //BAC
  int bac_nhits=0;
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

      //Tree
      if( Tu >0 || Td>0 ){
	event.bachitpat[i]= seg;
	bac_nhits++; 
      }
      event.bacua[i] = Au;
      event.bacut[i] = Tu;
      event.bacda[i] = Ad;
      event.bacdt[i] = Td;

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
      //       if( Td>0 ){
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
  event.bacnhits =bac_nhits;

  //TGT
  int tgt_nhits=0;
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

      //Tree
      if( Tu >0 || Td>0 ){
	event.tgthitpat[i]= seg;
	tgt_nhits++; 
      }
      event.tgta[i] = Au;
      event.tgtt[i] = Tu;

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
      //       if( Td>0 ){
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
  event.tgtnhits =tgt_nhits;

  //TOF
  int tof_nhits=0;
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

      //Tree
      if( Tu >0 || Td>0 ){
	event.tofhitpat[i]= seg;
	tof_nhits++; 
      }
      event.tofua[i] = Au;
      event.tofut[i] = Tu;
      event.tofda[i] = Ad;
      event.tofdt[i] = Td;

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
      if( Td>0 ){
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
  event.tofnhits =tof_nhits;

  //AC
  int ac1_nhits=0;
  int ac2_nhits=0;
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
	  //Tree
	  if( T >0 ){
	    event.ac1hitpat[i]= seg;
	    ac1_nhits++;
	  }
	  event.ac1a[i] = A;
	  event.ac1t[i] = T;

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
	  //Tree
	  if( T>0 ){
	    event.ac2hitpat[i]= seg+20;
	    ac2_nhits++; 
	  }
	  event.ac2a[i] = A;
	  event.ac2t[i] = T;

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
  event.ac1nhits =ac1_nhits;
  event.ac2nhits =ac2_nhits;

  //LC
  int lc_nhits=0;
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

      //Tree
      if( Tu >0 || Td>0 ){
	 event.lchitpat[i]= seg;
	 lc_nhits++; 
      } 
      event.lcua[i] = Au;
      event.lcut[i] = Tu;
      event.lcda[i] = Ad;
      event.lcdt[i] = Td;

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
      if( Td>0 ){
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
  event.lcnhits =lc_nhits;

  tree->Fill();
  return true;
}

void EventHodoMonitor::InitializeEvent( void )
{
  event.gcnhits  = -1;
  event.bh1nhits = -1;
  event.bh2nhits = -1;
  event.bacnhits = -1;
  event.tgtnhits = -1;
  event.tofnhits = -1;
  event.ac1nhits = -1;
  event.ac2nhits = -1;
  event.lcnhits  = -1;

  for( int it=0; it<MaxHits; it++){
    event.gchitpat[it]  = -1;
    event.bh1hitpat[it] = -1;
    event.bh2hitpat[it] = -1;
    event.bachitpat[it] = -1;
    event.tgthitpat[it] = -1;
    event.tofhitpat[it] = -1;
    event.ac1hitpat[it] = -1;
    event.ac2hitpat[it] = -1;
    event.lchitpat[it]  = -1;

    event.trigpat[it]  = -1;
  }

  for( int it=0; it<NumOfSegGC; it++){
    event.gca[it] = -999.0;
    event.gct[it] = -999.0;
  }

  for( int it=0; it<NumOfSegBH1; it++){
    event.bh1ua[it] = -999.0;
    event.bh1ut[it] = -999.0;
    event.bh1da[it] = -999.0;
    event.bh1dt[it] = -999.0;
  }

  for( int it=0; it<NumOfSegBH2; it++){
    event.bh2ua[it] = -999.0;
    event.bh2ut[it] = -999.0;
    event.bh2da[it] = -999.0;
    event.bh2dt[it] = -999.0;
  }

  for( int it=0; it<NumOfSegBAC; it++){
    event.bacua[it] = -999.0;
    event.bacut[it] = -999.0;
    event.bacda[it] = -999.0;
    event.bacdt[it] = -999.0;
  }

  for( int it=0; it<NumOfSegTGT; it++){
    event.tgta[it] = -999.0;
    event.tgtt[it] = -999.0;
  }

  for( int it=0; it<NumOfSegTOF; it++){
    event.tofua[it] = -999.0;
    event.tofut[it] = -999.0;
    event.tofda[it] = -999.0;
    event.tofdt[it] = -999.0;
  }

  for( int it=0; it<NumOfSegAC; it++){
    event.ac1a[it] = -999.0;
    event.ac1t[it] = -999.0;
    event.ac2a[it] = -999.0;
    event.ac2t[it] = -999.0;
  }

  for( int it=0; it<NumOfSegLC; it++){
    event.lcua[it] = -999.0;
    event.lcut[it] = -999.0;
    event.lcda[it] = -999.0;
    event.lcdt[it] = -999.0;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }
}

bool EventHodoMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventHodoMonitor;
}

const int NbinAdc = 1024;
const double MinAdc  =    0.;
const double MaxAdc  = 4096.;

const int NbinTdc = 1024;
const double MinTdc  =    0.;
const double MaxTdc  = 4096.;

bool ConfMan:: InitializeHistograms()
{
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
    title1 << "AC-1" << i << " Adc";
    HB1( ACHid+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title2 << "AC-2" << i << " Adc";
    HB1( ACHid+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title3 << "AC-1" << i << " Tdc";
    HB1( ACHid+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
    title4 << "AC-2" << i << " Tdc";
    HB1( ACHid+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );

    std::stringstream title5, title6;
    title5 << "AC-1" << i << " Adc(w Tdc)";
    HB1( ACHid+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title6 << "AC-2" << i << " Adc(w Tdc)";
    HB1( ACHid+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );

    std::stringstream title7, title8;
    title7 << "AC-1" << i << " Adc(w/o Tdc)";
    HB1( ACHid+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    title8 << "AC-2" << i << " Adc(w/o Tdc)";
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

 //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //GC
  tree->Branch("gcnhits",   &event.gcnhits,   "gcnhits/I");
  tree->Branch("gchitpat",   event.gchitpat,  "gchitpat[64]/I");
  tree->Branch("gca",        event.gca,       "gca[1]/D");
  tree->Branch("gct",        event.gct,       "gct[1]/D");

  //BH1
  tree->Branch("bh1nhits",   &event.bh1nhits,   "bh1nhits/I");
  tree->Branch("bh1hitpat",   event.bh1hitpat,  "bh1hitpat[64]/I");
  tree->Branch("bh1ua",       event.bh1ua,      "bh1ua[11]/D");
  tree->Branch("bh1ut",       event.bh1ut,      "bh1ut[11]/D");
  tree->Branch("bh1da",       event.bh1da,      "bh1da[11]/D");
  tree->Branch("bh1dt",       event.bh1dt,      "bh1dt[11]/D");

  //BH2
  tree->Branch("bh2nhits",   &event.bh2nhits,   "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,  "bh2hitpat[64]/I");
  tree->Branch("bh2ua",       event.bh2ua,      "bh2ua[7]/D");
  tree->Branch("bh2ut",       event.bh2ut,      "bh2ut[7]/D");
  tree->Branch("bh2da",       event.bh2da,      "bh2da[7]/D");
  tree->Branch("bh2dt",       event.bh2dt,      "bh2dt[7]/D");

  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  "bachitpat[64]/I");
  tree->Branch("bacua",       event.bacua,      "bacua[10]/D");
  tree->Branch("bacut",       event.bacut,      "bacut[10]/D");
  tree->Branch("bacda",       event.bacda,      "bacda[10]/D");
  tree->Branch("bacdt",       event.bacdt,      "bacdt[10]/D");

  //TGT
  tree->Branch("tgtnhits",   &event.tgtnhits,   "tgtnhits/I");
  tree->Branch("tgthitpat",   event.tgthitpat,  "tgthitpat[64]/I");
  tree->Branch("tgta",        event.tgta,       "tgta[3]/D");
  tree->Branch("tgtt",        event.tgtt,       "tgtt[3]/D");

  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  "tofhitpat[64]/I");
  tree->Branch("tofua",       event.tofua,      "tofua[32]/D");
  tree->Branch("tofut",       event.tofut,      "tofut[32]/D");
  tree->Branch("tofda",       event.tofda,      "tofda[32]/D");
  tree->Branch("tofdt",       event.tofdt,      "tofdt[32]/D");

  //AC
  tree->Branch("ac1nhits",   &event.ac1nhits,   "ac1nhits/I");
  tree->Branch("ac2nhits",   &event.ac2nhits,   "ac2nhits/I");
  tree->Branch("ac1hitpat",   event.ac1hitpat,  "ac1hitpat[64]/I");
  tree->Branch("ac2hitpat",   event.ac2hitpat,  "ac2hitpat[64]/I");
  tree->Branch("ac1a",        event.ac1a,       "ac1a[20]/D");
  tree->Branch("ac1t",        event.ac1t,       "ac1t[20]/D");
  tree->Branch("ac2a",        event.ac2a,       "ac2a[20]/D");
  tree->Branch("ac2t",        event.ac2t,       "ac2t[20]/D");

  //LC
  tree->Branch("lcnhits",   &event.lcnhits,   "lcnhits/I");
  tree->Branch("lchitpat",   event.lchitpat,  "lchitpat[64]/I");
  tree->Branch("lcua",       event.lcua,      "lcua[28]/D");
  tree->Branch("lcut",       event.lcut,      "lcut[28]/D");
  tree->Branch("lcda",       event.lcda,      "lcda[28]/D");
  tree->Branch("lcdt",       event.lcdt,      "lcdt[28]/D");

  //Misc
  tree->Branch("trigpat",    event.trigpat,   "trigpat[64]/I");
  tree->Branch("trigflag",   event.trigflag,   "trigflag[10]/I");

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
