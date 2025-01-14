/*
  UserHodoscopeSP0.cc
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

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "HodoCluster.hh"
#include "BH2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "BH2Cluster.hh"

static const int GCHid =   9000;
static const int BH1Hid = 10000;
static const int BH2Hid = 20000;
static const int BACHid = 30000;
static const int TGTHid = 40000;
static const int TOFHid = 50000;
static const int ACHid  = 60000;
static const int LCHid  = 70000;
static const int SP0Hid = 80000;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventHodoscopeSP0 
  : public VEvent
{

private:
  RawData *rawData;
  HodoAnalyzer *hodoAna;

public:
  EventHodoscopeSP0();
  ~EventHodoscopeSP0();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

EventHodoscopeSP0::EventHodoscopeSP0()
  : VEvent(),
    rawData(0),
    hodoAna(new HodoAnalyzer)
{
}

EventHodoscopeSP0::~EventHodoscopeSP0()
{
  delete hodoAna;
  if (rawData) delete rawData;
}

bool EventHodoscopeSP0::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

#ifndef MaxHits 
#define MaxHits 32
#endif
#ifndef MaxHits2 
#define MaxHits2 20
#endif

//For Tree
struct Event{
  int trigpat[NumOfMisc];
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

  double bh1mt[NumOfSegBH1];
  double bh1de[NumOfSegBH1];
  double bh2mt[NumOfSegBH2];
  double bh2de[NumOfSegBH2];
  double tofmt[NumOfSegTOF];
  double tofde[NumOfSegTOF];
  double lcmt[NumOfSegLC];
  double lcde[NumOfSegLC];

  int sp0nhits[NumOfLayersSP0];
  int sp0hitpat[NumOfLayersSP0][MaxHits];
  double sp0ua[NumOfLayersSP0][NumOfSegSP0];
  double sp0ut[NumOfLayersSP0][NumOfSegSP0];
  double sp0da[NumOfLayersSP0][NumOfSegSP0];
  double sp0dt[NumOfLayersSP0][NumOfSegSP0];

//   int nhitclb;
//   double btof[MaxHits2];

//   int nhitcls;
//   double stof[MaxHits2];
};
static Event event;

const double AdcMin  = 0.;
const double AdcMax  = 4096.;
const double TdcMin  = 0.;
const double TdcMax  = 4096.;

bool EventHodoscopeSP0::ProcessingNormal()
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
      event.trigflag[seg-1] = T;
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
	event.gchitpat[gc_nhits] = seg;
	gc_nhits++; 
      }
      event.gca[seg-1] = A;
      event.gct[seg-1] = T;

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
	event.bh1hitpat[bh1_nhits] = seg;
	bh1_nhits++; 
      }
      event.bh1ua[seg-1] = Au;
      event.bh1ut[seg-1] = Tu;
      event.bh1da[seg-1] = Ad;
      event.bh1dt[seg-1] = Td;

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
	event.bh2hitpat[bh2_nhits]= seg;
	bh2_nhits++; 
      }
      event.bh2ua[seg-1] = Au;
      event.bh2ut[seg-1] = Tu;
      event.bh2da[seg-1] = Ad;
      event.bh2dt[seg-1] = Td;

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
	event.bachitpat[bac_nhits]= seg;
	bac_nhits++; 
      }
      event.bacua[seg-1] = Au;
      event.bacut[seg-1] = Tu;
      event.bacda[seg-1] = Ad;
      event.bacdt[seg-1] = Td;

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
	event.tgthitpat[tgt_nhits]= seg;
	tgt_nhits++; 
      }
      event.tgta[seg-1] = Au;
      event.tgtt[seg-1] = Tu;

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
	event.tofhitpat[tof_nhits]= seg;
	tof_nhits++; 
      }
      event.tofua[seg-1] = Au;
      event.tofut[seg-1] = Tu;
      event.tofda[seg-1] = Ad;
      event.tofdt[seg-1] = Td;

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
    for( int layer=1; layer<=NumOfLayersAc; ++layer ){
      const HodoRHitContainer &cont=rawData->GetACRawHC(layer);
      int nh=cont.size();
      if( layer==1 ) HF1( ACHid, double(nh) );
      if( layer==2 ) HF1( ACHid, double(nh)+20. );
      int nh1=0, nh2=0;
      for( int i=0; i<nh; ++i ){
	HodoRawHit *hit=cont[i];
	int seg=hit->SegmentId()+1;
	std::string name="AC";
	HF1( name.c_str(), seg-0.5 );
	int A=hit->GetAdc1();
	int T=hit->GetTdc1();

	if( layer==1 ){
	  //Tree
	  if( T>0 ){
	    event.ac1hitpat[ac1_nhits]= seg;
	    ac1_nhits++;
	  }
	  event.ac1a[seg-1] = A;
	  event.ac1t[seg-1] = T;

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
	if( layer==2 ){
	  //Tree
	  if( T>0 ){
	    event.ac2hitpat[ac2_nhits]= seg;
	    ac2_nhits++; 
	  }
	  event.ac2a[seg-1] = A;
	  event.ac2t[seg-1] = T;

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
      if( layer==1 ){
	HF1( ACHid+2, double(nh1) ); 
	HF1( ACHid+4, double(nh2) );
      }
      if( layer==2 ){
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
	 event.lchitpat[lc_nhits]= seg;
	 lc_nhits++; 
      } 
      event.lcua[seg-1] = Au;
      event.lcut[seg-1] = Tu;
      event.lcda[seg-1] = Ad;
      event.lcdt[seg-1] = Td;

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

  //SP0
  int sp0_nhits=0; 
  {  
    for( int layer=1; layer<=NumOfLayersSP0; ++layer ){
      const HodoRHitContainer &cont=rawData->GetSP0RawHC(layer);
      int nh=cont.size();
      HF1( SP0Hid+1000*layer, double(nh) );
      int nh1=0, nh2=0;
      event.sp0nhits[layer-1] = nh;
      for( int i=0; i<nh; ++i ){
	HodoRawHit *hit=cont[i];
	int seg=hit->SegmentId()+1;
	std::string name="SP0";
	HF1( name.c_str()+layer, seg-0.5 );
	int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
	int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

	//Tree
	if( Tu >0 || Td>0 ){
	  event.sp0hitpat[layer-1][sp0_nhits]= seg;
	  sp0_nhits++; 
	} 
	event.sp0ua[layer-1][seg-1] = Au;
	event.sp0ut[layer-1][seg-1] = Tu;
	event.sp0da[layer-1][seg-1] = Ad;
	event.sp0dt[layer-1][seg-1] = Td;
	
	//Up
	HF1( SP0Hid+1000*layer+100*seg+1, double(Au) );
	if( Tu>0 ){
	  HF1( SP0Hid+1000*layer+100*seg+3, double(Tu) );
	  HF1( SP0Hid+1000*layer+100*seg+5, double(Au) );
	}
	else{
	  HF1( SP0Hid+1000*layer+100*seg+7, double(Au) );
	}
	//Down
	HF1( SP0Hid+1000*layer+100*seg+2, double(Ad) );
	if( Td>0 ){
	  HF1( SP0Hid+1000*layer+100*seg+4, double(Td) );
	  HF1( SP0Hid+1000*layer+100*seg+6, double(Ad) );
	}
	else{
	  HF1( SP0Hid+1000*layer+100*seg+7, double(Ad) );
	}
	
	//HitPat
	if( Tu>0 || Td>0 ){
	  ++nh1; HF1( SP0Hid+1000*layer+3, seg-0.5 );
	}
	if( Tu>0 && Td>0 ){
	  ++nh2; HF1( SP0Hid+1000*layer+5, seg-0.5 );
	}
      }
      HF1( SP0Hid+1000*layer+2, double(nh1) ); 
      HF1( SP0Hid+1000*layer+4, double(nh2) );
    }
  }

  //**************************************************************************
  //******************NormalizedData

  //BH1
  hodoAna->DecodeBH1Hits( rawData );
  {
    int nh=hodoAna->GetNHitsBH1();
    HF1( BH1Hid+10, double(nh) );
    int nh2=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitBH1(i);
      event.bh1mt[i] = hit->MeanTime();
      event.bh1de[i] = hit->DeltaE();
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( BH1Hid+11, seg-0.5 );
      double au=hit->GetAUp(), ad=hit->GetADown();
      double tu=hit->GetTUp(), td=hit->GetTDown();
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double mt=hit->MeanTime(), cmt=hit->CMeanTime();
      double de=hit->DeltaE();
      HF1( BH1Hid+100*seg+11, tu ); HF1( BH1Hid+100*seg+12, td );
      HF1( BH1Hid+100*seg+13, mt ); HF1( BH1Hid+100*seg+14, au );
      HF1( BH1Hid+100*seg+15, ad ); HF1( BH1Hid+100*seg+16, de );
      HF1( BH1Hid+100*seg+17, ctu ); HF1( BH1Hid+100*seg+18, ctd );
      HF1( BH1Hid+100*seg+19, cmt ); HF1( BH1Hid+100*seg+20, ctu-ctd );
      HF2( BH1Hid+100*seg+21, tu, au );  HF2( BH1Hid+100*seg+22, td, ad );
      HF2( BH1Hid+100*seg+23, ctu, au ); HF2( BH1Hid+100*seg+24, ctd, ad );
      HF1( BH1Hid+12, cmt ); 
      HF1( BH1Hid+13, de ); 
      if( de>0.5 ){
	++nh2; HF1( BH1Hid+15, seg-0.5 ); 
	HF1( BH1Hid+16, cmt );
      }
    }
    HF1( BH1Hid+14, double(nh2) );
    for( int i1=0; i1<nh; ++i1 ){
      Hodo2Hit *hit1=hodoAna->GetHitBH1(i1);
      if(!hit1 || hit1->DeltaE()<=0.5 ) continue;
      int seg1=hit1->SegmentId()+1;
      for( int i2=0; i2<nh; ++i2 ){
	if( i1==i2 ) continue;
	Hodo2Hit *hit2=hodoAna->GetHitBH1(i2);
	if( !hit2 || hit2->DeltaE()<=0.5 ) continue;
	int seg2=hit2->SegmentId()+1;
	double ct1=hit1->CMeanTime(), ct2=hit2->CMeanTime();
	HF2( BH1Hid+21, seg1-0.5, seg2-0.5 );
	HF2( BH1Hid+22, ct1, ct2 );
	HF1( BH1Hid+23, ct2-ct1 );
	if( fabs(ct2-ct1)<2.0 ){
	  HF2( BH1Hid+24, seg1-0.5, seg2-0.5 );
	}
      }
    }
  }

  // BH2
  hodoAna->DecodeBH2Hits( rawData );
  {
    int nh=hodoAna->GetNHitsBH2();
    HF1( BH2Hid+10, double(nh) );
    int nh2=0;
    for( int i=0; i<nh; ++i ){
      BH2Hit *hit=hodoAna->GetHitBH2(i);
      event.bh2mt[i] = hit->MeanTime();
      event.bh2de[i] = hit->DeltaE();
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( BH2Hid+11, seg-0.5 );
      double au=hit->GetAUp(), ad=hit->GetADown();
      double tu=hit->GetTUp(), td=hit->GetTDown();
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double mt=hit->MeanTime(), cmt=hit->CMeanTime();
      double de=hit->DeltaE();
      HF1( BH2Hid+100*seg+11, tu ); HF1( BH2Hid+100*seg+12, td );
      HF1( BH2Hid+100*seg+13, mt ); HF1( BH2Hid+100*seg+14, au );
      HF1( BH2Hid+100*seg+15, ad ); HF1( BH2Hid+100*seg+16, de );
      HF1( BH2Hid+100*seg+17, ctu ); HF1( BH2Hid+100*seg+18, ctd );
      HF1( BH2Hid+100*seg+19, cmt ); HF1( BH2Hid+100*seg+20, ctu-ctd );
      HF2( BH2Hid+100*seg+21, tu, au );  HF2( BH2Hid+100*seg+22, td, ad );
      HF2( BH2Hid+100*seg+23, ctu, au ); HF2( BH2Hid+100*seg+24, ctd, ad );
      HF1( BH2Hid+12, cmt ); HF1( BH2Hid+13, de ); 
      if( de>0.5 ){
	++nh2; HF1( BH2Hid+15, seg-0.5 ); HF1( BH2Hid+16, cmt );
      }
    }
    HF1( BH2Hid+14, double(nh2) );
    for( int i1=0; i1<nh; ++i1 ){
      Hodo2Hit *hit1=hodoAna->GetHitBH2(i1);
      if(!hit1 || hit1->DeltaE()<=0.5 ) continue;
      int seg1=hit1->SegmentId()+1;
      for( int i2=0; i2<nh; ++i2 ){
	if( i1==i2) continue;
	Hodo2Hit *hit2=hodoAna->GetHitBH2(i2);
	if(!hit2 || hit2->DeltaE()<=0.5 ) continue;
	int seg2=hit2->SegmentId()+1;
	double ct1=hit1->CMeanTime(), ct2=hit2->CMeanTime();
	HF2( BH2Hid+21, seg1-0.5, seg2-0.5 );
	HF2( BH2Hid+22, ct1, ct2 );
	HF1( BH2Hid+23, ct2-ct1 );
	if( fabs(ct2-ct1)<2.0 ){
	  HF2( BH2Hid+24, seg1-0.5, seg2-0.5 );
	}
      }
    }

    int nc=hodoAna->GetNClustersBH2();
    HF1( BH2Hid+30, double(nc) );
    int nc2=0;
    for( int i=0; i<nc; ++i ){
      BH2Cluster *cluster=hodoAna->GetClusterBH2(i);
      if(!cluster) continue;
      int cs=cluster->ClusterSize();
      double ms=cluster->MeanSeg()+1,
	cmt=cluster->CMeanTime(), de=cluster->DeltaE();
      HF1( BH2Hid+31, double(cs) );
      HF1( BH2Hid+32, ms-0.5 );
      HF1( BH2Hid+33, cmt ); HF1( BH2Hid+34, de );
      if( de>0.5 ){
	++nc2; HF1( BH2Hid+36, cmt );
      }

      for( int i2=0; i2<nc; ++i2 ){
	if( i2==i ) continue;
	BH2Cluster *cl2=hodoAna->GetClusterBH2(i2);
	if(!cl2) continue;
	double ms2=cl2->MeanSeg()+1, cmt2=cl2->CMeanTime(),
	  de2=cl2->DeltaE();
	if( de<=0.5 || de2<=0.5 ) continue;
	HF2( BH2Hid+41, ms-0.5, ms2-0.5 );
	HF2( BH2Hid+42, cmt, cmt2 );
	HF1( BH2Hid+43, cmt2-cmt );
	if( fabs(cmt2-cmt)<2.0 ){
	  HF2( BH2Hid+44, ms-0.5, ms2-0.5 );
	}
      }
    }
    HF1( BH2Hid+35, double(nc2) );
  }

  // BH1 with BH2 gate
  {
    int nhbh2=hodoAna->GetNHitsBH2();
    //    if( nhbh2==1 ){
    if( nhbh2 ){
      int seg2=hodoAna->GetHitBH2(0)->SegmentId()+1;
      double mt2=hodoAna->GetHitBH2(0)->CTime0();
      int nh=hodoAna->GetNHitsBH1();
      for( int i=0; i<nh; ++i ){
	Hodo2Hit *hit=hodoAna->GetHitBH1(i);
	if(!hit) continue;
	int seg1=hit->SegmentId()+1;
	double tu1=hit->GetTUp(), td1=hit->GetTDown();
	double mt1=hit->MeanTime();
	HF1( BH1Hid+100*seg1+1100+21+seg2*10, tu1 );
	HF1( BH1Hid+100*seg1+1100+22+seg2*10, td1 );
	HF1( BH1Hid+100*seg1+1100+23+seg2*10, mt1 );

	//For BH1vsBH2 Correlation
	HF2( BH1Hid+100*seg1+2200+21+seg2*10, tu1, mt2 );
	HF2( BH1Hid+100*seg1+2200+22+seg2*10, td1, mt2 );
	HF2( BH1Hid+100*seg1+2200+23+seg2*10, mt1, mt2 );
      }
    }
    for( int i2=0; i2<nhbh2; ++i2 ){
      int seg2=hodoAna->GetHitBH2(i2)->SegmentId()+1;
      double mt0=hodoAna->GetHitBH2(i2)->CTime0();
      int nhbh1=hodoAna->GetNHitsBH1();
      for( int i=0; i<nhbh1; ++i ){
	Hodo2Hit *hit=hodoAna->GetHitBH1(i);
	if(!hit) continue;
	int seg1=hit->SegmentId()+1;
	double mt1=hit->MeanTime();
	HF1( BH1Hid+100*seg1+1100+24+seg2*10, mt1-mt0 );
	HF1( BH1Hid+100*seg1+2200+104, mt1-mt0 );
      }
    }
  }

  // BH1-BH2
  {
    int ncbh1=hodoAna->GetNClustersBH1();
    int ncbh2=hodoAna->GetNClustersBH2();
    int nhit_clb=0;
    for( int i1=0; i1<ncbh1; ++i1 ){
      HodoCluster *clbh1=hodoAna->GetClusterBH1(i1);
      if(!clbh1) continue;
      double seg1=clbh1->MeanSeg()+1, mt1=clbh1->CMeanTime();
      double de1=clbh1->DeltaE();
      for(int i2=0; i2<ncbh2; ++i2 ){
	BH2Cluster *clbh2=hodoAna->GetClusterBH2(i2);
	if(!clbh2) continue;
	double seg2=clbh2->MeanSeg()+1, t0=clbh2->CTime0();
	double de2=clbh2->DeltaE();
	HF1( 201, mt1-t0 );
	HF2( 202, seg1-0.5, seg2-0.5 );
	//For BH1vsBH2 Correlation
	HF2( 203, t0, mt1 );
	HF1( 204, mt1 );
	HF1( 205, t0 );

// 	event.btof[nhit_clb]  = mt1-t0;
// 	nhit_clb++;

	if( de1>0.5 && de2>0.5 ){
	  HF1( 211, mt1-t0 );
	  HF2( 212, seg1-0.5, seg2-0.5 );
	}
      }
    }
//     event.nhitclb = nhit_clb;
  }

  // BH1-BH2 PHC
  {
    int nh1=hodoAna->GetNHitsBH1();
    int nh2=hodoAna->GetNHitsBH2();
    for( int i2=0; i2<nh2; ++i2 ){
      BH2Hit *hit2=hodoAna->GetHitBH2(i2);
      int seg2=hit2->SegmentId()+1;
      double au2=hit2->GetAUp(), ad2=hit2->GetADown();
      double tu2=hit2->GetTUp(), td2=hit2->GetTDown();
      double ctu2=hit2->GetCTUp(), ctd2=hit2->GetCTDown();
      double time0=hit2->Time0(), ctime0=hit2->CTime0();
      double tofs=ctime0-0.5*(ctu2+ctd2);
      for( int i1=0; i1<nh1; ++i1 ){
	Hodo2Hit *hit1=hodoAna->GetHitBH1(i1);
	int seg1=hit1->SegmentId()+1;
	double au1=hit1->GetAUp(), ad1=hit1->GetADown();
	double tu1=hit1->GetTUp(), td1=hit1->GetTDown();
	double ctu1=hit1->GetCTUp(), ctd1=hit1->GetCTDown();
	double cmt1=hit1->CMeanTime();
	
	HF2( 100*seg1+BH1Hid+81, 2.*ctime0-ctu1-ctd1, au1 );
	HF2( 100*seg1+BH1Hid+82, 2.*ctime0-ctu1-ctd1, ad1 );
	HF2( 100*seg1+BH1Hid+83, 2.*ctime0-tu1-ctd1, au1 );
	HF2( 100*seg1+BH1Hid+84, 2.*ctime0-ctu1-td1, ad1 );
	
	HF2( 100*seg2+BH2Hid+81, 2.*(cmt1-tofs)-ctu2-ctd2, au2 );
	HF2( 100*seg2+BH2Hid+82, 2.*(cmt1-tofs)-ctu2-ctd2, ad2 );
	HF2( 100*seg2+BH2Hid+83, 2.*(cmt1-tofs)-tu2-ctd2, au2 );
	HF2( 100*seg2+BH2Hid+84, 2.*(cmt1-tofs)-ctu2-td2, ad2 );
      }
    }
  }

  // TOF
  hodoAna->DecodeTOFHits( rawData );
  {
    int nh=hodoAna->GetNHitsTOF();
    HF1( TOFHid+10, double(nh) );
    int nh2=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitTOF(i);
      event.tofmt[i] =hit->MeanTime();
      event.tofde[i] =hit->DeltaE();
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( TOFHid+11, seg-0.5 );
      double au=hit->GetAUp(), ad=hit->GetADown();
      double tu=hit->GetTUp(), td=hit->GetTDown();
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double mt=hit->MeanTime(), cmt=hit->CMeanTime();
      double de=hit->DeltaE();
      HF1( TOFHid+100*seg+11, tu ); HF1( TOFHid+100*seg+12, td );
      HF1( TOFHid+100*seg+13, mt ); HF1( TOFHid+100*seg+14, au );
      HF1( TOFHid+100*seg+15, ad ); HF1( TOFHid+100*seg+16, de );
      HF1( TOFHid+100*seg+17, ctu ); HF1( TOFHid+100*seg+18, ctd );
      HF1( TOFHid+100*seg+19, cmt ); HF1( TOFHid+100*seg+20, ctu-ctd );
      HF2( TOFHid+100*seg+21, tu, au );  HF2( TOFHid+100*seg+22, td, ad );
      HF2( TOFHid+100*seg+23, ctu, au ); HF2( TOFHid+100*seg+24, ctd, ad );
      HF1( TOFHid+12, cmt ); HF1( TOFHid+13, de );
      if( de>0.5 ){
	HF1( TOFHid+15, seg-0.5 );
	++nh2;
      }
      HF1( TOFHid+14, double(nh2) );
    }

    for( int i1=0; i1<nh; ++i1 ){
      Hodo2Hit *hit1=hodoAna->GetHitTOF(i1);
      if(!hit1 || hit1->DeltaE()<=0.5 ) continue;
      int seg1=hit1->SegmentId()+1;
      for( int i2=0; i2<nh; ++i2 ){
	if( i1==i2) continue;
	Hodo2Hit *hit2=hodoAna->GetHitTOF(i2);
	if(!hit2 || hit2->DeltaE()<=0.5 ) continue;
	int seg2=hit2->SegmentId()+1;
	double ct1=hit1->CMeanTime(), ct2=hit2->CMeanTime();
	HF2( TOFHid+21, seg1-0.5, seg2-0.5 );
	HF2( TOFHid+22, ct1, ct2 );
	HF1( TOFHid+23, ct2-ct1 );
	if( fabs(ct2-ct1)<3.0 ){
	  HF2( TOFHid+24, seg1-0.5, seg2-0.5 );
	}
      }
    }

    int nc=hodoAna->GetNClustersTOF();
    HF1( TOFHid+30, double(nc) );
    for( int i=0; i<nc; ++i ){
      HodoCluster *cluster=hodoAna->GetClusterTOF(i);
      if(!cluster) continue;
      int cs=cluster->ClusterSize();
      double ms=cluster->MeanSeg()+1,
	cmt=cluster->CMeanTime(), de=cluster->DeltaE();
      HF1( TOFHid+31, double(cs) );
      HF1( TOFHid+32, ms-0.5 );
      HF1( TOFHid+33, cmt ); HF1( TOFHid+34, de );
    }
  }

  // LC
  hodoAna->DecodeLCHits( rawData );
  {
    int nh=hodoAna->GetNHitsLC();
    HF1( LCHid+10, double(nh) );
    double nh2=0;
    for( int i=0; i<nh; ++i ){
      Hodo2Hit *hit=hodoAna->GetHitLC(i);
      event.lcmt[i] =hit->MeanTime();
      event.lcde[i] =hit->DeltaE();
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( LCHid+11, seg-0.5 );
      double au=hit->GetAUp(), ad=hit->GetADown();
      double tu=hit->GetTUp(), td=hit->GetTDown();
      double ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      double mt=hit->MeanTime(), cmt=hit->CMeanTime();
      double de=hit->DeltaE();
      HF1( LCHid+100*seg+11, tu ); HF1( LCHid+100*seg+12, td );
      HF1( LCHid+100*seg+13, mt ); HF1( LCHid+100*seg+14, au );
      HF1( LCHid+100*seg+15, ad ); HF1( LCHid+100*seg+16, de );
      HF1( LCHid+100*seg+17, ctu ); HF1( LCHid+100*seg+18, ctd );
      HF1( LCHid+100*seg+19, cmt ); HF1( LCHid+100*seg+20, ctu-ctd );
      HF2( LCHid+100*seg+21, tu, au );  HF2( LCHid+100*seg+22, td, ad );
      HF2( LCHid+100*seg+23, ctu, au ); HF2( LCHid+100*seg+24, ctd, ad );
      HF1( LCHid+12, cmt ); HF1( LCHid+13, de );
      if( de>0.2 ){
	HF1( LCHid+15, seg-0.5 );
	++nh2;
      }
      HF1( LCHid+14, double(nh2) );
    }
    for( int i1=0; i1<nh; ++i1 ){
      Hodo2Hit *hit1=hodoAna->GetHitLC(i1);
      if(!hit1 || hit1->DeltaE()<=0.2 ) continue;
      int seg1=hit1->SegmentId()+1;
      for( int i2=0; i2<nh; ++i2 ){
	if( i1==i2 ) continue;
	Hodo2Hit *hit2=hodoAna->GetHitLC(i2);
	if(!hit2 || hit2->DeltaE()<=0.2 ) continue;
	int seg2=hit2->SegmentId()+1;
	double ct1=hit1->CMeanTime(), ct2=hit2->CMeanTime();
	HF2( LCHid+21, seg1-0.5, seg2-0.5 );
	HF2( LCHid+22, ct1, ct2 );
	HF1( LCHid+23, ct2-ct1 );
	if( fabs(ct2-ct1)<5.0 ){
	  HF2( LCHid+24, seg1-0.5, seg2-0.5 );
	}
      }
    }

    int nc=hodoAna->GetNClustersLC();
    HF1( LCHid+30, double(nc) );
    for( int i=0; i<nc; ++i ){
      HodoCluster *cluster=hodoAna->GetClusterLC(i);
      if(!cluster) continue;
      int cs=cluster->ClusterSize();
      double ms=cluster->MeanSeg()+1,
	cmt=cluster->CMeanTime(), de=cluster->DeltaE();
      HF1( LCHid+31, double(cs) );
      HF1( LCHid+32, ms-0.5 );
      HF1( LCHid+33, cmt ); HF1( LCHid+34, de );
    }
  }

  //TOF-LC
  {
    int nhTof=hodoAna->GetNHitsTOF();
    int nhLc =hodoAna->GetNHitsLC();
    int nhit_cls=0;
    for( int iTof=0; iTof<nhTof; ++iTof ){
      for( int iLc=0; iLc<nhLc; ++iLc ){
	Hodo2Hit *hitTof=hodoAna->GetHitTOF(iTof);
	Hodo2Hit *hitLc =hodoAna->GetHitLC(iLc);
	if( !hitTof || !hitLc ) continue;
	int segTof=hitTof->SegmentId()+1, segLc=hitLc->SegmentId()+1;
	double cmtTof=hitTof->CMeanTime(), cmtLc=hitLc->CMeanTime();
	double deTof=hitTof->DeltaE(), deLc=hitLc->DeltaE();
	HF2( 101, segTof-0.5, segLc-0.5 );
	HF1( 103, cmtLc-cmtTof );
	HF2( 105, deTof, deLc );
	HF2( 106, cmtLc-cmtTof, deTof );
	HF2( 107, cmtLc-cmtTof, deLc );
// 	event.stof[nhit_cls] = cmtLc-cmtTof;
// 	nhit_cls++;

	if( deTof>0.6 && deLc>0.2 ){
	  HF2( 102, segTof-0.5, segLc-0.5 );
	  HF1( 104, cmtLc-cmtTof );
	}
      }
    }
//     event.nhitcls = nhit_cls;
  }

  //AC
  hodoAna->DecodeACHits( rawData );
  for( int layer=1; layer<=2; ++layer ){
    int nh=hodoAna->GetNHitsAC(layer);
    HF1( ACHid+10+50*(layer-1), double(nh) );
    for( int i=0; i<nh; ++i ){
      Hodo1Hit *hit=hodoAna->GetHitAC(layer,i);
      if(!hit) continue;
      int seg=hit->SegmentId()+1;
      HF1( ACHid+11+50*(layer-1), seg-0.5 );
      double a=hit->GetA(), t=hit->GetT(), ct=hit->GetCT();
      HF1( ACHid+100*seg+11+50*(layer-1), t);  
      HF1( ACHid+100*seg+12+50*(layer-1), a); 
      HF1( ACHid+100*seg+13+50*(layer-1), ct); 
    }
  }

  tree->Fill();
  return true;
}

void EventHodoscopeSP0::InitializeEvent( void )
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

    event.bh1mt[it]  = -999.0;
    event.bh1de[it]  = -999.0;
  }

  for( int it=0; it<NumOfSegBH2; it++){
    event.bh2ua[it] = -999.0;
    event.bh2ut[it] = -999.0;
    event.bh2da[it] = -999.0;
    event.bh2dt[it] = -999.0;

    event.bh2mt[it]  = -999.0;
    event.bh2de[it]  = -999.0;
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

    event.tofmt[it]  = -999.0;
    event.tofde[it]  = -999.0;
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

    event.lcmt[it] = -999.0;
    event.lcde[it] = -999.0;
  }

  for( int it=0; it<NumOfMisc; it++){
    event.trigpat[it] = -1;
    event.trigflag[it] = -1;
  }


  //SP0
  for( int it=0; it<NumOfLayersSP0; it++){
    event.sp0nhits[it] = -1;
    for( int jt=0; jt<NumOfSegSP0; jt++){
      event.sp0ua[it][jt] = -999.0;
      event.sp0ut[it][jt] = -999.0;
      event.sp0da[it][jt] = -999.0;
      event.sp0dt[it][jt] = -999.0;
    }
    for( int jt=0; jt<MaxHits; jt++){
      event.sp0hitpat[it][jt] = -1;
    }
  }

//   event.nhitclb  = -1;
//   event.nhitcls = -1;

//   for( int it=0; it<MaxHits2; it++){
//     event.btof[it]   = -999.0;

//     event.stof[it]   = -999.0;
//   }
}

bool EventHodoscopeSP0::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventHodoscopeSP0;
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
  // Rawdata
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

  //BH1 Normalized
  HB1( BH1Hid+10, "#Hits BH1[Hodo]",  NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+11, "Hitpat BH1[Hodo]", NumOfSegBH1,   0., double(NumOfSegBH1)   );
  HB1( BH1Hid+12, "CMeanTime BH1", 200, -10., 10. );
  HB1( BH1Hid+13, "dE BH1", 200, -0.5, 4.5 );
  HB1( BH1Hid+14, "#Hits BH1[HodoGood]",  NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+15, "Hitpat BH1[HodoGood]", NumOfSegBH1,   0., double(NumOfSegBH1)   );
  HB1( BH1Hid+16, "CMeanTime BH1[HodoGood]", 200, -10., 10. );

  for( int i=1; i<=NumOfSegBH1; ++i ){
    std::stringstream title1, title2, title3;
    title1 << "BH1-" << i << "Up Time";
    HB1( BH1Hid+100*i+11, title1.str().c_str(), 200, -10., 10. );
    title2 << "BH1-" << i << "Down Time";
    HB1( BH1Hid+100*i+12, title2.str().c_str(), 200, -10., 10. );
    title3 << "BH1-" << i << "MeanTime";
    HB1( BH1Hid+100*i+13, title3.str().c_str(), 200, -10., 10. );

    std::stringstream title4, title5, title6;
    title4 << "BH1-" << i << "Up dE";
    HB1( BH1Hid+100*i+14, title4.str().c_str(), 200, -0.5, 4.5 );
    title5 << "BH1-" << i << "Down dE";
    HB1( BH1Hid+100*i+15, title5.str().c_str(), 200, -0.5, 4.5 );
    title6 << "BH1-" << i << "dE";
    HB1( BH1Hid+100*i+16, title6.str().c_str(), 200, -0.5, 4.5 );

    std::stringstream title7, title8, title9;
    title7 << "BH1-" << i << "Up CTime";
    HB1( BH1Hid+100*i+17, title7.str().c_str(), 200, -10., 10. );
    title8 << "BH1-" << i << "Down CTime";
    HB1( BH1Hid+100*i+18, title8.str().c_str(), 200, -10., 10. );
    title9 << "BH1-" << i << "CMeanTime";
    HB1( BH1Hid+100*i+19, title9.str().c_str(), 200, -10., 10. );

    std::stringstream title10;
    title10 << "BH1-" << i << "Tup-Tdown";
    HB1( BH1Hid+100*i+20, title10.str().c_str(), 200, -5.0, 5.0 );

    std::stringstream title11, title12, title13, title14;
    title11 << "BH1-" << i << "Up dE%Time"; 
    HB2( BH1Hid+100*i+21, title11.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title12 << "BH1-" << i << "Down dE%Time"; 
    HB2( BH1Hid+100*i+22, title12.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title13 << "BH1-" << i << "Up dE%CTime"; 
    HB2( BH1Hid+100*i+23, title13.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title14 << "BH1-" << i << "Down dE%CTime"; 
    HB2( BH1Hid+100*i+24, title14.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );

    //For BH1vsBH2 Correlation
    std::stringstream title15, title16, title17;
    title15 << "BH1-" << i << "BH2 Up MT%MT"; 
    HB2( BH1Hid+100*i+2200+21, title15.str().c_str(), 100, -10., 10., 100, -10., 10. );
    title16 << "BH1-" << i << "BH2 Down MT%MT"; 
    HB2( BH1Hid+100*i+2200+22, title16.str().c_str(), 100, -10., 10., 100, -10., 10. );
    title17 << "BH1-" << i << "BH2 MeanTime MT%MT"; 
    HB2( BH1Hid+100*i+2200+23, title17.str().c_str(), 100, -10., 10., 100, -10., 10. );
  }
  for( int i=1; i<=NumOfSegBH1; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BH1-" << i << "Up Time [BH2]";
    HB1( BH1Hid+1100+100*i+31, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+41, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+51, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+61, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+71, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+81, title1.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+91, title1.str().c_str(), 200, -10., 10. );
    title2 << "BH1-" << i << "Down Time [BH2]";
    HB1( BH1Hid+1100+100*i+32, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+42, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+52, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+62, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+72, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+82, title2.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+92, title2.str().c_str(), 200, -10., 10. );
    title3 << "BH1-" << i << "MeanTime [BH2]";
    HB1( BH1Hid+1100+100*i+33, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+43, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+53, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+63, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+73, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+83, title3.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+93, title3.str().c_str(), 200, -10., 10. );
    title4 << "BH1-" << i << "MeanTime-BH2MeamTime";
    HB1( BH1Hid+1100+100*i+34, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+44, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+54, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+64, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+74, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+84, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+1100+100*i+94, title4.str().c_str(), 200, -10., 10. );
    HB1( BH1Hid+100*i+2200+104, title4.str().c_str(), 200, -10., 10. );
  }

  HB2( BH1Hid+21, "BH1HitPat%BH1HitPat[HodoGood]", NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH1,   0., double(NumOfSegBH1) ); 
  HB2( BH1Hid+22, "CMeanTimeBH1%CMeanTimeBH1[HodoGood]", 
       100, -5., 5., 100, -5., 5. );
  HB1( BH1Hid+23, "TDiff BH1[HodoGood]", 200, -10., 10. ); 
  HB2( BH1Hid+24, "BH1HitPat%BH1HitPat[HodoGood2]", NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH1,   0., double(NumOfSegBH1) ); 

  HB1( BH1Hid+30, "#Clusters BH1", NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+31, "ClusterSize BH1", 5, 0., 5. ); 
  HB1( BH1Hid+32, "HitPat Cluster BH1", 2*NumOfSegBH1, 0., double(NumOfSegBH1) );
  HB1( BH1Hid+33, "CMeamTime Cluster BH1", 200, -10., 10. ); 
  HB1( BH1Hid+34, "DeltaE Cluster BH1", 100, -0.5, 4.5 );
  HB1( BH1Hid+35, "#Clusters BH1(AdcGood)", NumOfSegBH1+1, 0., double(NumOfSegBH1+1) );
  HB1( BH1Hid+36, "CMeamTime Cluster BH1(AdcGood)", 200, -10., 10. ); 

  HB2( BH1Hid+41, "BH1ClP%BH1ClP",  NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH1,   0., double(NumOfSegBH1) ); 
  HB2( BH1Hid+42, "CMeanTimeBH1%CMeanTimeBH1[Cluster]",
       100, -5., 5., 100, -5., 5. );
  HB1( BH1Hid+43, "TDiff BH1[Cluster]", 200, -10., 10. ); 
  HB2( BH1Hid+44, "BH1ClP%BH1ClP[AdcGood]",  NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH1,   0., double(NumOfSegBH1) ); 

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

  HB1( BH2Hid+10, "#Hits BH2[Hodo]",  NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+11, "Hitpat BH2[Hodo]", NumOfSegBH2,   0., double(NumOfSegBH2)   );
  HB1( BH2Hid+12, "CMeanTime BH2", 200, -10., 10. );
  HB1( BH2Hid+13, "dE BH2", 200, -0.5, 4.5 );
  HB1( BH2Hid+14, "#Hits BH2[HodoGood]",  NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+15, "Hitpat BH2[HodoGood]", NumOfSegBH2,   0., double(NumOfSegBH2)   );
  HB1( BH2Hid+16, "CMeanTime BH2[HodoGood]", 200, -10., 10. );

  for( int i=1; i<=NumOfSegBH2; ++i ){
    std::stringstream title1, title2, title3;
    title1 << "BH2-" << i << "Up Time";
    HB1( BH2Hid+100*i+11, title1.str().c_str(), 200, -10., 10. );
    title2 << "BH2-" << i << "Down Time";
    HB1( BH2Hid+100*i+12, title2.str().c_str(), 200, -10., 10. );
    title3 << "BH2-" << i << "MeanTime";
    HB1( BH2Hid+100*i+13, title3.str().c_str(), 200, -10., 10. );

    std::stringstream title4, title5, title6;
    title4 << "BH2-" << i << "Up dE";
    HB1( BH2Hid+100*i+14, title4.str().c_str(), 200, -0.5, 4.5 );
    title5 << "BH2-" << i << "Down dE";
    HB1( BH2Hid+100*i+15, title5.str().c_str(), 200, -0.5, 4.5 );
    title6 << "BH2-" << i << "dE";
    HB1( BH2Hid+100*i+16, title6.str().c_str(), 200, -0.5, 4.5 );

    std::stringstream title7, title8, title9;
    title1 << "BH2-" << i << "Up CTime";
    HB1( BH2Hid+100*i+17, title7.str().c_str(), 200, -10., 10. );
    title2 << "BH2-" << i << "Down CTime";
    HB1( BH2Hid+100*i+18, title8.str().c_str(), 200, -10., 10. );
    title3 << "BH2-" << i << "CMeanTime";
    HB1( BH2Hid+100*i+19, title9.str().c_str(), 200, -10., 10. );

    std::stringstream title10;
    title10 << "BH2-" << i << "Tup-Tdown";
    HB1( BH2Hid+100*i+20, title10.str().c_str(), 200, -5.0, 5.0 );

    std::stringstream title11, title12, title13, title14;
    title11 << "BH2-" << i << "Up dE%Time"; 
    HB2( BH2Hid+100*i+21, title11.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title12 << "BH2-" << i << "Down dE%Time"; 
    HB2( BH2Hid+100*i+22, title12.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title13 << "BH2-" << i << "Up dE%CTime"; 
    HB2( BH2Hid+100*i+23, title13.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
    title14 << "BH2-" << i << "Down dE%CTime"; 
    HB2( BH2Hid+100*i+24, title14.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 );
  }

  HB2( BH2Hid+21, "BH2HitPat%BH2HitPat[HodoGood]", NumOfSegBH2,   0., double(NumOfSegBH2),
       NumOfSegBH2,   0., double(NumOfSegBH2) ); 
  HB2( BH2Hid+22, "CMeanTimeBH2%CMeanTimeBH2[HodoGood]",
       100, -2.5, 2.5, 100, -2.5, 2.5 );
  HB1( BH2Hid+23, "TDiff BH2[HodoGood]", 200, -10., 10. ); 
  HB2( BH2Hid+24, "BH2HitPat%BH2HitPat[HodoGood2]", NumOfSegBH2,   0., double(NumOfSegBH2),
       NumOfSegBH2,   0., double(NumOfSegBH2) ); 

  HB1( BH2Hid+30, "#Clusters BH2", NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+31, "ClusterSize BH2", 5, 0., 5. ); 
  HB1( BH2Hid+32, "HitPat Cluster BH2", 2*NumOfSegBH2, 0., double(NumOfSegBH2) );
  HB1( BH2Hid+33, "CMeamTime Cluster BH2", 200, -10., 10. ); 
  HB1( BH2Hid+34, "DeltaE Cluster BH2", 100, -0.5, 4.5 );
  HB1( BH2Hid+35, "#Clusters BH2(ADCGood)", NumOfSegBH2+1, 0., double(NumOfSegBH2+1) );
  HB1( BH2Hid+36, "CMeamTime Cluster BH2(ADCGood)", 200, -10., 10. ); 

  HB2( BH2Hid+41, "BH2ClP%BH2ClP", NumOfSegBH2,   0., double(NumOfSegBH2),
       NumOfSegBH2,   0., double(NumOfSegBH2) ); 
  HB2( BH2Hid+42, "CMeanTimeBH2%CMeanTimeBH2[Cluster]",
       100, -2.5, 2.5, 100, -2.5, 2.5 );
  HB1( BH2Hid+43, "TDiff BH2[Cluster]", 200, -10., 10. ); 
  HB2( BH2Hid+44, "BH2ClP%BH2ClP(ADCGood)", NumOfSegBH2,   0., double(NumOfSegBH2),
       NumOfSegBH2,   0., double(NumOfSegBH2) ); 


  HB1( 201, "TimeDif BH1-BH2", 200, -10., 10. );
  HB2( 202, "SegBH2%SegBH1",NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH2,   0., double(NumOfSegBH2) );
  //For BH1vsBH2 Corr
  HB2( 203, "MTBH2%MTBH1", 200, -10., 10., 200, -10., 10. );
  HB1( 204, "MTBH1", 200, -10., 10.);
  HB1( 205, "MTBH2", 200, -10., 10.);

  HB1( 211, "TimeDif BH1-BH2(GoodAdc)", 200, -10., 10. );
  HB2( 212, "SegBH2%SegBH1(GoodAdc)",NumOfSegBH1,   0., double(NumOfSegBH1),
       NumOfSegBH2,   0., double(NumOfSegBH2) );
  

  // BH1-BH2 PHC
  for( int i=1; i<=NumOfSegBH1; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BH1-" << i << " U dE%CT-TOF";
    title2 << "BH1-" << i << " D dE%CT-TOF";
    title3 << "BH1-" << i << " U dE%T-TOF";
    title4 << "BH1-" << i << " D dE%T-TOF";
    HB2( BH1Hid+100*i+81, title1.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH1Hid+100*i+82, title2.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH1Hid+100*i+83, title3.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH1Hid+100*i+84, title4.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
  }  
  for( int i=1; i<=NumOfSegBH2; ++i ){
    std::stringstream title1, title2, title3, title4;
    title1 << "BH2-" << i << " U dE%CT-TOF";
    title2 << "BH2-" << i << " D dE%CT-TOF";
    title3 << "BH2-" << i << " U dE%T-TOF";
    title4 << "BH2-" << i << " D dE%T-TOF";
    HB2( BH2Hid+100*i+81, title1.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH2Hid+100*i+82, title2.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH2Hid+100*i+83, title3.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
    HB2( BH2Hid+100*i+84, title4.str().c_str(), 100, -10., 10., 100, -0.5, 4.5 ); 
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

  HB1( TOFHid+10, "#Hits Tof[Hodo]",  NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+11, "Hitpat Tof[Hodo]", NumOfSegTOF,   0., double(NumOfSegTOF)   );
  HB1( TOFHid+12, "CMeanTime Tof", 500, -5., 45. );
  HB1( TOFHid+13, "dE Tof", 200, -0.5, 4.5 );
  HB1( TOFHid+14, "#Hits Tof[HodoGood]",  NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+15, "Hitpat Tof[HodoGood]", NumOfSegTOF,   0., double(NumOfSegTOF)   );

  for( int i=1; i<=NumOfSegTOF; ++i ){
    std::stringstream title1, title2, title3;
    title1 << "TOF-" << i << "Up Time";
    HB1( TOFHid+100*i+11, title1.str().c_str(), 500, -5., 45. );
    title2 << "TOF-" << i << "Down Time";
    HB1( TOFHid+100*i+12, title2.str().c_str(), 500, -5., 45. );
    title3 << "TOF-" << i << "MeanTime";
    HB1( TOFHid+100*i+13, title3.str().c_str(), 500, -5., 45. );

    std::stringstream title4, title5, title6;
    title4 << "TOF-" << i << "Up dE";
    HB1( TOFHid+100*i+14, title4.str().c_str(), 200, -0.5, 4.5 );
    title5 << "TOF-" << i << "Down dE";
    HB1( TOFHid+100*i+15, title5.str().c_str(), 200, -0.5, 4.5 );
    title6 << "TOF-" << i << "dE";
    HB1( TOFHid+100*i+16, title6.str().c_str(), 200, -0.5, 4.5 );

    std::stringstream title7, title8, title9;
    title7 << "TOF-" << i << "Up CTime";
    HB1( TOFHid+100*i+17, title7.str().c_str(), 500, -5., 45. );
    title8 << "TOF-" << i << "Down CTime";
    HB1( TOFHid+100*i+18, title8.str().c_str(), 500, -5., 45. );
    title9 << "TOF-" << i << "CMeanTime";
    HB1( TOFHid+100*i+19, title9.str().c_str(), 500, -5., 45. );

    std::stringstream title10;
    title10 << "TOF-" << i << "Tup-Tdown";
    HB1( TOFHid+100*i+20, title10.str().c_str(), 200, -10.0, 10.0 );
  }

  HB2( TOFHid+21, "TofHitPat%TofHitPat[HodoGood]", NumOfSegTOF,   0., double(NumOfSegTOF),
       NumOfSegTOF,   0., double(NumOfSegTOF) ); 
  HB2( TOFHid+22, "CMeanTimeTof%CMeanTimeTof[HodoGood]",
       120, 10., 40., 120, 10., 40. );
  HB1( TOFHid+23, "TDiff Tof[HodoGood]", 200, -10., 10. ); 
  HB2( TOFHid+24, "TofHitPat%TofHitPat[HodoGood2]", NumOfSegTOF,   0., double(NumOfSegTOF),
       NumOfSegTOF,   0., double(NumOfSegTOF) ); 

  HB1( TOFHid+30, "#Clusters Tof", NumOfSegTOF+1, 0., double(NumOfSegTOF+1) );
  HB1( TOFHid+31, "ClusterSize Tof", 5, 0., 5. ); 
  HB1( TOFHid+32, "HitPat Cluster Tof", 2*NumOfSegTOF, 0., double(NumOfSegTOF) );
  HB1( TOFHid+33, "CMeamTime Cluster Tof", 500, -5., 45. ); 
  HB1( TOFHid+34, "DeltaE Cluster Tof", 100, -0.5, 4.5 );

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

  HB1( ACHid+10, "#Hits AC1[Hodo]",     NumOfSegAC+1, 0., double(NumOfSegAC+1) );
  HB1( ACHid+11, "Hitpat AC1[Hodo]",    NumOfSegAC,   0., double(NumOfSegAC)   );
  HB1( ACHid+60, "#Hits AC2[Hodo]",     NumOfSegAC+1, 0., double(NumOfSegAC+1) );
  HB1( ACHid+61, "Hitpat AC2[Hodo]",    NumOfSegAC,   0., double(NumOfSegAC)   );

  for( int i=1; i<=NumOfSegAC; ++i ){
    std::stringstream title1, title2;
    title1 << "Ac-1-" << i << " Time";
    HB1( ACHid+100*i+11, title1.str().c_str(), 500, -5., 45. );
    title2 << "Ac-2-" << i << " Time";
    HB1( ACHid+100*i+61, title2.str().c_str(), 500, -5., 45. );
    std::stringstream title3, title4;
    title3 << "Ac-1-" << i << " dE";
    HB1( ACHid+100*i+12, title3.str().c_str(), 200, -0.5, 4.5 ); 
    title4 << "Ac-2-" << i << " dE";
    HB1( ACHid+100*i+62, title4.str().c_str(), 200, -0.5, 4.5 ); 
    std::stringstream title5, title6;
    title5 << "Ac-1-" << i << " CTime";
    HB1( ACHid+100*i+13, title5.str().c_str(), 500, -5., 45. );
    title6 << "Ac-2-" << i << " CTime";
    HB1( ACHid+100*i+63, title6.str().c_str(), 500, -5., 45. );
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

  HB1( LCHid+10, "#Hits Lc[Hodo]",  NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+11, "Hitpat Lc[Hodo]", NumOfSegLC,   0., double(NumOfSegLC)   );
  HB1( LCHid+12, "CMeanTime Lc", 500, -5., 45. );
  HB1( LCHid+13, "dE Lc", 200, -0.5, 4.5 );
  HB1( LCHid+14, "#Hits Lc[HodoGood]",  NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+15, "Hitpat Lc[HodoGood]", NumOfSegLC,   0., double(NumOfSegLC)   );

  for( int i=1; i<=NumOfSegLC; ++i ){
    std::stringstream title1, title2, title3;
    title1 << "LC-" << i << "Up Time";
    HB1( LCHid+100*i+11, title1.str().c_str(), 500, -5., 45. );
    title2 << "LC-" << i << "Down Time";
    HB1( LCHid+100*i+12, title2.str().c_str(), 500, -5., 45. );
    title3 << "LC-" << i << "MeanTime";
    HB1( LCHid+100*i+13, title3.str().c_str(), 500, -5., 45. );

    std::stringstream title4, title5, title6;
    title4 << "LC-" << i << "Up dE";
    HB1( LCHid+100*i+14, title4.str().c_str(), 200, -0.5, 4.5 );
    title5 << "LC-" << i << "Down dE";
    HB1( LCHid+100*i+15, title5.str().c_str(), 200, -0.5, 4.5 );
    title6 << "LC-" << i << "dE";
    HB1( LCHid+100*i+16, title6.str().c_str(), 200, -0.5, 4.5 );

    std::stringstream title7, title8, title9;
    title7 << "LC-" << i << "Up CTime";
    HB1( LCHid+100*i+17, title7.str().c_str(), 500, -5., 45. );
    title8 << "LC-" << i << "Down CTime";
    HB1( LCHid+100*i+18, title8.str().c_str(), 500, -5., 45. );
    title9 << "LC-" << i << "CMeanTime";
    HB1( LCHid+100*i+19, title9.str().c_str(), 500, -5., 45. );

    std::stringstream title10;
    title10 << "LC-" << i << "Tup-Tdown";
    HB1( LCHid+100*i+20, title10.str().c_str(), 300, -15.0, 15.0 );
  }
 
  HB2( LCHid+21, "LcHitPat%LcHitPat[HodoGood]", NumOfSegLC,   0., double(NumOfSegLC),
       NumOfSegLC,   0., double(NumOfSegLC) ); 
  HB2( LCHid+22, "CMeanTimeLc%CMeanTimeLc[HodoGood]",
       120, 10., 40., 120, 10., 40. );
  HB1( LCHid+23, "TDiff Lc[HodoGood]", 200, -10., 10. ); 
  HB2( LCHid+24, "LcHitPat%LcHitPat[HodoGood2]", NumOfSegLC,   0., double(NumOfSegLC),
       NumOfSegLC,   0., double(NumOfSegLC) ); 

  HB1( LCHid+30, "#Clusters Lc", NumOfSegLC+1, 0., double(NumOfSegLC+1) );
  HB1( LCHid+31, "ClusterSize Lc", 5, 0., 5. ); 
  HB1( LCHid+32, "HitPat Cluster Lc", 2*NumOfSegLC, 0., double(NumOfSegLC) );
  HB1( LCHid+33, "CMeamTime Cluster Lc", 500, -5., 45. ); 
  HB1( LCHid+34, "DeltaE Cluster Lc", 100, -0.5, 4.5 );

  //TOF-LC
  HB2( 101, "SegLc%SegTof", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., NumOfSegLC );
  HB2( 102, "SegLc%SegTof[Good]", NumOfSegTOF, 0., double(NumOfSegTOF),
       NumOfSegLC, 0., NumOfSegLC );
  HB1( 103, "LcTime-TofTime", 300, -5., 25. );
  HB1( 104, "LcTime-TofTime[Good]", 300, -5., 25. );

  HB2( 105, "dELc%dETof", 100, -0.5, 4.5, 100, -0.5, 4.5 ); 
  HB2( 106, "dETof%LcTime-TofTime",100, -5., 20., 100, -0.5, 4.5 );
  HB2( 107, "dELc%LcTime-TofTime",100, -5., 20., 100, -0.5, 4.5 );

  // SP0
  for( int j=1; j<=NumOfLayersSP0; ++j ){
    HB1( SP0Hid+1000*j, "#Hits SP0",        NumOfSegSP0+1, 0., double(NumOfSegSP0+1) );
    HB1( SP0Hid+1000*j+1, "Hitpat SP0",       NumOfSegSP0,   0., double(NumOfSegSP0)   );
    HB1( SP0Hid+1000*j+2, "#Hits SP0(Tor)",   NumOfSegSP0+1, 0., double(NumOfSegSP0+1) );
    HB1( SP0Hid+1000*j+3, "Hitpat SP0(Tor)",  NumOfSegSP0,   0., double(NumOfSegSP0)   );
    HB1( SP0Hid+1000*j+4, "#Hits SP0(Tand)",  NumOfSegSP0+1, 0., double(NumOfSegSP0+1) );
    HB1( SP0Hid+1000*j+5, "Hitpat SP0(Tand)", NumOfSegSP0,   0., double(NumOfSegSP0)   );

    for( int i=1; i<=NumOfSegSP0; ++i ){
      std::stringstream title1, title2, title3, title4;
      title1 << "SP0-" << j << "-" << i << " UpAdc";
      HB1( SP0Hid+1000*j+100*i+1, title1.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title2 << "SP0-" << j << "-" << i << " DownAdc";
      HB1( SP0Hid+1000*j+100*i+2, title2.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title3 << "SP0-" << j << "-" << i << " UpTdc";
      HB1( SP0Hid+1000*j+100*i+3, title3.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      title4 << "SP0-" << j << "-" << i << " DownTdc";
      HB1( SP0Hid+1000*j+100*i+4, title4.str().c_str(), NbinTdc, MinTdc, MaxTdc );
      
      std::stringstream title5, title6;
      title5 << "SP0-" << j << "-" << i << " UpAdc(w Tdc)";
      HB1( SP0Hid+1000*j+100*i+5, title5.str().c_str(), NbinAdc, MinAdc, MaxAdc );
      title6 << "SP0-" << j << "-" << i << " DownAdc(w Tdc)";
      HB1( SP0Hid+1000*j+100*i+6, title6.str().c_str(), NbinAdc, MinAdc, MaxAdc );
      
      std::stringstream title7, title8;
      title7 << "SP0-" << j << "-" << i << " UpAdc(w/o Tdc)";
      HB1( SP0Hid+1000*j+100*i+7, title7.str().c_str(), NbinAdc, MinAdc, MaxAdc );
      title8 << "SP0-" << j << "-" << i << " DownAdc(w/o Tdc)";
      HB1( SP0Hid+1000*j+100*i+8, title8.str().c_str(), NbinAdc, MinAdc, MaxAdc );
    }
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //GC
  tree->Branch("gcnhits",   &event.gcnhits,   "gcnhits/I");
  tree->Branch("gchitpat",   event.gchitpat,  "gchitpat[32]/I");
  tree->Branch("gca",        event.gca,       "gca[1]/D");
  tree->Branch("gct",        event.gct,       "gct[1]/D");

  //BH1
  tree->Branch("bh1nhits",   &event.bh1nhits,   "bh1nhits/I");
  tree->Branch("bh1hitpat",   event.bh1hitpat,  "bh1hitpat[32]/I");
  tree->Branch("bh1ua",       event.bh1ua,      "bh1ua[11]/D");
  tree->Branch("bh1ut",       event.bh1ut,      "bh1ut[11]/D");
  tree->Branch("bh1da",       event.bh1da,      "bh1da[11]/D");
  tree->Branch("bh1dt",       event.bh1dt,      "bh1dt[11]/D");

  //BH2
  tree->Branch("bh2nhits",   &event.bh2nhits,   "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,  "bh2hitpat[32]/I");
  tree->Branch("bh2ua",       event.bh2ua,      "bh2ua[7]/D");
  tree->Branch("bh2ut",       event.bh2ut,      "bh2ut[7]/D");
  tree->Branch("bh2da",       event.bh2da,      "bh2da[7]/D");
  tree->Branch("bh2dt",       event.bh2dt,      "bh2dt[7]/D");

  //BAC
  tree->Branch("bacnhits",   &event.bacnhits,   "bacnhits/I");
  tree->Branch("bachitpat",   event.bachitpat,  "bachitpat[32]/I");
  tree->Branch("bacua",       event.bacua,      "bacua[10]/D");
  tree->Branch("bacut",       event.bacut,      "bacut[10]/D");
  tree->Branch("bacda",       event.bacda,      "bacda[10]/D");
  tree->Branch("bacdt",       event.bacdt,      "bacdt[10]/D");

  //TGT
  tree->Branch("tgtnhits",   &event.tgtnhits,   "tgtnhits/I");
  tree->Branch("tgthitpat",   event.tgthitpat,  "tgthitpat[32]/I");
  tree->Branch("tgta",        event.tgta,       "tgta[4]/D");
  tree->Branch("tgtt",        event.tgtt,       "tgtt[4]/D");

  //TOF
  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  "tofhitpat[32]/I");
  tree->Branch("tofua",       event.tofua,      "tofua[32]/D");
  tree->Branch("tofut",       event.tofut,      "tofut[32]/D");
  tree->Branch("tofda",       event.tofda,      "tofda[32]/D");
  tree->Branch("tofdt",       event.tofdt,      "tofdt[32]/D");

  //AC
  tree->Branch("ac1nhits",   &event.ac1nhits,   "ac1nhits/I");
  tree->Branch("ac2nhits",   &event.ac2nhits,   "ac2nhits/I");
  tree->Branch("ac1hitpat",   event.ac1hitpat,  "ac1hitpat[32]/I");
  tree->Branch("ac2hitpat",   event.ac2hitpat,  "ac2hitpat[32]/I");
  tree->Branch("ac1a",        event.ac1a,       "ac1a[20]/D");
  tree->Branch("ac1t",        event.ac1t,       "ac1t[20]/D");
  tree->Branch("ac2a",        event.ac2a,       "ac2a[20]/D");
  tree->Branch("ac2t",        event.ac2t,       "ac2t[20]/D");

  //LC
  tree->Branch("lcnhits",   &event.lcnhits,   "lcnhits/I");
  tree->Branch("lchitpat",   event.lchitpat,  "lchitpat[32]/I");
  tree->Branch("lcua",       event.lcua,      "lcua[28]/D");
  tree->Branch("lcut",       event.lcut,      "lcut[28]/D");
  tree->Branch("lcda",       event.lcda,      "lcda[28]/D");
  tree->Branch("lcdt",       event.lcdt,      "lcdt[28]/D");

  //Misc
  tree->Branch("trigpat",    event.trigpat,   "trigpat[20]/I");
  tree->Branch("trigflag",   event.trigflag,  "trigflag[20]/I");

  //Normalized data
  tree->Branch("bh1mt",     event.bh1mt,     "bh1mt[11]/D");  
  tree->Branch("bh1de",     event.bh1de,     "bh1de[11]/D");  
  tree->Branch("bh2mt",     event.bh2mt,     "bh2mt[7]/D");
  tree->Branch("bh2de",     event.bh2de,     "bh2de[7]/D");  
  tree->Branch("tofmt",     event.tofmt,     "tofmt[32]/D");  
  tree->Branch("tofde",     event.tofde,     "tofde[32]/D");  
  tree->Branch("lcmt",      event.lcmt,      "lcmt[28]/D");
  tree->Branch("lcde",      event.lcde,     " lcde[28]/D");  

  //SP0
  tree->Branch("sp0nhits",   &event.sp0nhits,   "sp0nhits[8]/I");
  tree->Branch("sp0hitpat",   event.sp0hitpat,  "sp0hitpat[8][5]/I");
  tree->Branch("sp0ua",       event.sp0ua,      "sp0ua[8][5]/D");
  tree->Branch("sp0ut",       event.sp0ut,      "sp0ut[8][5]/D");
  tree->Branch("sp0da",       event.sp0da,      "sp0da[8][5]/D");
  tree->Branch("sp0dt",       event.sp0dt,      "sp0dt[8][5]/D");

//   tree->Branch("nhitclb",  &event.nhitclb,   "nhitclb/I");
//   tree->Branch("btof",      event.btof,      "btof[nhitclb]/D");  

//   tree->Branch("nhitcls",  &event.nhitcls,   "nhitcls/I");
//   tree->Branch("stof",      event.stof,      "stof[nhitclb]/D");  

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
//   if( CMapFileName_!="" )
//     CMapManager_ = new CMapMan(CMapFileName_);
//   if(CMapManager_) CMapManager_->Initialize();

  if( HodoParamFileName_!="" )
    HodoParamManager_ = new HodoParamMan(HodoParamFileName_);
  if(HodoParamManager_) HodoParamManager_->Initialize();

  if( HodoPHCFileName_!="" )
    HodoPHCManager_ = new HodoPHCMan(HodoPHCFileName_);
  if( HodoPHCManager_ ) HodoPHCManager_->Initialize();

//   if( ScalerDefinitionFileName_!="" )
//     ScalerAnalyzer_ = new ScalerAna(ScalerDefinitionFileName_);
//   else
//     ScalerAnalyzer_ = new ScalerAna();
//   if(ScalerAnalyzer_) ScalerAnalyzer_->Initialize(); 

  return true;
}
