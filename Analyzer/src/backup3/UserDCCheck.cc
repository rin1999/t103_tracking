/*
  UserDCCheck.cc
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
#include "DCRawHit.hh"

#include "HodoRawHit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "BH2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "BH2Cluster.hh"

const double TdcLow  =    0.;
const double TdcHigh = 1000.;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventDCCheck 
  : public VEvent
{

private:
  RawData *rawData;
  HodoAnalyzer *hodoAna;

public:
  EventDCCheck();
  ~EventDCCheck();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventDCCheck::EventDCCheck()
  : VEvent(),
    rawData(0),
    hodoAna(new HodoAnalyzer)
{
}

EventDCCheck::~EventDCCheck()
{
  delete hodoAna;
  if (rawData) delete rawData;
}

bool EventDCCheck::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

#ifndef MaxHits 
#define MaxHits 100
#endif
#ifndef MaxHits2 
#define MaxHits2 64
#endif
#ifndef MaxHits3 
#define MaxHits3 1000
#endif

struct Event{
  int trigtype;
  int trigflag[NumOfMisc];

  int bc1nhits[NumOfLayersBc];
  int bc1hitpat[NumOfLayersBc][MaxHits];

  int bc2nhits[NumOfLayersBc];
  int bc2hitpat[NumOfLayersBc][MaxHits];

  int bc3nhits[NumOfLayersBc];
  int bc3hitpat[NumOfLayersBc][MaxHits];

  int bc4nhits[NumOfLayersBc];
  int bc4hitpat[NumOfLayersBc][MaxHits];

  int sdc1nhits[NumOfLayersSdc];
  int sdc1hitpat[NumOfLayersSdc][MaxHits];

  int sdc2nhits[NumOfLayersSdc];
  int sdc2hitpat[NumOfLayersSdc][MaxHits];

  int sdc3nhits[NumOfLayersSdc];
  int sdc3hitpat[NumOfLayersSdc][MaxHits];

  int sdc4nhits[NumOfLayersSdc];
  int sdc4hitpat[NumOfLayersSdc][MaxHits];

  int bh1nhits;
  int bh1hitpat[MaxHits3];

  int bh2nhits;
  int bh2hitpat[MaxHits3];

  int tofnhits;
  int tofhitpat[MaxHits2];

  int nhitclb;
  double bh1mt[MaxHits3];
  double bh1de[MaxHits3];
  double bh2mt[MaxHits3];
  double bh2de[MaxHits3];
  double btof[MaxHits3];

  int nhitcls;
  double tofmt[MaxHits3];
  double tofde[MaxHits3];
  double lcmt[MaxHits3];
  double lcde[MaxHits3];
  double stof[MaxHits3];
};
static Event event;

bool EventDCCheck::ProcessingNormal()
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


  //BH1
  int bh1_nhits=0;
  {
    const HodoRHitContainer &cont=rawData->GetBH1RawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="BH1";
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Tree
      if( Tu>0 || Td>0 ){
	event.bh1hitpat[i] = seg;
	bh1_nhits++; 
      }
    }
  }
  event.bh1nhits = bh1_nhits;

  //BH2
  int bh2_nhits=0;
  {
    const HodoRHitContainer &cont=rawData->GetBH2RawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="BH2";
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Tree
      if( Tu>0 || Td>0 ){
	event.bh2hitpat[i] = seg;
	bh2_nhits++; 
      }
    }
  }
  event.bh2nhits = bh2_nhits;

  //TOF
  int tof_nhits=0;
  {
    const HodoRHitContainer &cont=rawData->GetTOFRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg=hit->SegmentId()+1;
      std::string name="TOF";
      int Au=hit->GetAdcUp(), Ad=hit->GetAdcDown();
      int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();

      //Tree
      if( Tu>0 || Td>0 ){
	event.tofhitpat[i] = seg;
	tof_nhits++; 
      }
    }
  }
  event.tofnhits = tof_nhits;

  //BH1-BH2 for PID
  hodoAna->DecodeBH1Hits( rawData );
  hodoAna->DecodeBH2Hits( rawData );
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

	event.bh1mt[nhit_clb] = mt1;
	event.bh1de[nhit_clb] = de1;
	event.bh2mt[nhit_clb] = t0;
	event.bh2de[nhit_clb] = de2;
	event.btof[nhit_clb]  = mt1-t0;
	nhit_clb++;
      }
    }
    event.nhitclb = nhit_clb;
  }

  //TOF
  hodoAna->DecodeTOFHits( rawData );
  {
    int nc=hodoAna->GetNClustersTOF();
    for( int i=0; i<nc; ++i ){
      HodoCluster *cluster=hodoAna->GetClusterTOF(i);
      if(!cluster) continue;
      int cs=cluster->ClusterSize();
      double ms=cluster->MeanSeg()+1,
	cmt=cluster->CMeanTime(), de=cluster->DeltaE();
      
      event.tofmt[nc] =cmt;
      event.tofde[nc] =de;
    }
  }
  
  //LC
  hodoAna->DecodeLCHits( rawData );
  {
    int nc=hodoAna->GetNClustersLC();
    for( int i=0; i<nc; ++i ){
      HodoCluster *cluster=hodoAna->GetClusterLC(i);
      if(!cluster) continue;
      int cs=cluster->ClusterSize();
      double ms=cluster->MeanSeg()+1,
	cmt=cluster->CMeanTime(), de=cluster->DeltaE();

      event.lcmt[nc] =cmt;
      event.lcde[nc] =de;
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

	event.stof[nhit_cls] = cmtLc-cmtTof;
	nhit_cls++;
      }
    }
    event.nhitcls = nhit_cls;
  }


  //**************************************************************************
  //******************RawData

  //BC1&BC2
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      const DCRHitContainer &contIn =rawData->GetBcInRawHC(layer);
      int nhIn=contIn.size();
      if( layer<=NumOfLayersBc ){
	HF1( 100*layer+0, nhIn );
	event.bc1nhits[layer-1] = nhIn;
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();
	  
	  HF1( 100*layer+1, wire-0.5 );
	  event.bc1hitpat[layer-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( 0< tdc && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*layer+3, 10*tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 100*layer+2, wire-0.5 );
	  }
	  HF1( 100*layer+4, 10*tdc1st );
	}
      }
      else{
	event.bc2nhits[layer-NumOfLayersBc-1] = nhIn;
	HF1( 100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 100*layer+1, wire-0.5 );	  
	  event.bc2hitpat[layer-NumOfLayersBc-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( 0< tdc && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*layer+3, 10*tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 100*layer+2, wire-0.5 );
	  }
	  HF1( 100*layer+4, 10*tdc1st );
	}

      }
    }
  }  

  //BC3&BC4
  int bcIn_planeNum=12;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      const DCRHitContainer &contOut =rawData->GetBcOutRawHC(layer);
      int nhOut=contOut.size();
      if( layer<=NumOfLayersBc ){
	event.bc3nhits[layer-1] = nhOut;
	HF1( 100*(layer+bcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+bcIn_planeNum)+1, wire-0.5 );
	  event.bc3hitpat[layer-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*(layer+bcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 100*(layer+bcIn_planeNum)+2, wire-0.5 );
	  }
	  HF1( 100*(layer+bcIn_planeNum)+4, tdc1st );
	}
      }
      else{
	event.bc4nhits[layer-NumOfLayersBc-1] = nhOut;
	HF1( 100*(layer+bcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+bcIn_planeNum)+1, wire-0.5 );	  
	  event.bc4hitpat[layer-NumOfLayersBc-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*(layer+bcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 100*(layer+bcIn_planeNum)+2, wire-0.5 );
	  }
	  HF1( 100*(layer+bcIn_planeNum)+4, tdc1st );
	}
      }
    }
  }  

  //SDC1&SDC2
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      const DCRHitContainer &contIn =rawData->GetSdcInRawHC(layer);
      int nhIn=contIn.size();
      if( layer<=NumOfLayersSdc-2 ){
	event.sdc1nhits[layer-1] = nhIn;
	HF1( 3000+100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();
	  
	  HF1( 3000+100*layer+1, wire-0.5 );
	  event.sdc1hitpat[layer-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 3000+100*layer+2, wire-0.5 );
	  }
	  HF1( 3000+100*layer+4, tdc1st );
	}
      }
      else{
	event.sdc2nhits[layer-(NumOfLayersSdc-2)-1] = nhIn;
	HF1( 3000+100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*layer+1, wire-0.5 );	  
	  event.sdc2hitpat[layer-(NumOfLayersSdc-2)-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 3000+100*layer+2, wire-0.5 );
	  }
	  HF1( 3000+100*layer+4, tdc1st );
	}
      }
    }
  }

  //SDC3&SDC4
  int sdcIn_planeNum=10;
  {
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      const DCRHitContainer &contOut =rawData->GetSdcOutRawHC(layer);
      int nhOut=contOut.size();
      if( layer<=NumOfLayersSdc ){
	event.sdc3nhits[layer-1] = nhOut;
	HF1( 3000+100*(layer+sdcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+sdcIn_planeNum)+1, wire-0.5 );
	  event.sdc3hitpat[layer-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( tdc ) tdc1st=tdc;
	    HF1( 3000+100*(layer+sdcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 3000+100*(layer+sdcIn_planeNum)+2, wire-0.5 );
	  }
	  HF1( 3000+100*(layer+sdcIn_planeNum)+4, tdc1st );
	}
      }
      else{
	event.sdc4nhits[layer-NumOfLayersSdc-1] = nhOut;
	HF1( 3000+100*(layer+sdcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+sdcIn_planeNum)+1, wire-0.5 );	  
	  event.sdc4hitpat[layer-NumOfLayersSdc-1][i] = wire;
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( tdc ) tdc1st=tdc;
	    HF1( 3000+100*(layer+sdcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ){
	    HF1( 3000+100*(layer+sdcIn_planeNum)+2, wire-0.5 );
	  }
	  HF1( 3000+100*(layer+sdcIn_planeNum)+4, tdc1st );
	}
      }
    }
  }
  
  tree->Fill();
  return true;
}

bool EventDCCheck::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}


void EventDCCheck::InitializeEvent( void )
{
  event.trigtype = -1;
  for( int it=0; it<NumOfMisc; it++){
    event.trigflag[it] = -1;
  }

  for( int it=0; it<NumOfLayersBc; it++){
    event.bc1nhits[it]    = -1;
    event.bc2nhits[it]    = -1;
    event.bc3nhits[it]    = -1;
    event.bc4nhits[it]    = -1;

    event.sdc1nhits[it]   = -1;
    event.sdc2nhits[it]   = -1;
    event.sdc3nhits[it]   = -1;
    event.sdc4nhits[it]   = -1;
  }

  for( int it=0; it<NumOfLayersBc; it++){
    for( int jt=0; jt<MaxHits; jt++){
      event.bc1hitpat[it][jt] = -1;
      event.bc2hitpat[it][jt] = -1;
      event.bc3hitpat[it][jt] = -1;
      event.bc4hitpat[it][jt] = -1;
      
      event.sdc1hitpat[it][jt] = -1;
      event.sdc2hitpat[it][jt] = -1;
      event.sdc3hitpat[it][jt] = -1;
      event.sdc4hitpat[it][jt] = -1;
    }
  }

  event.bh1nhits = -1;
  event.bh2nhits = -1;
  event.tofnhits = -1;

  event.nhitclb  = -1;
  event.nhitcls  = -1;

  for( int it=0; it<MaxHits3; it++){
    event.bh1hitpat[it] = -1;
    event.bh2hitpat[it] = -1;
 
    event.bh1mt[it]  = -999.0;
    event.bh1de[it]  = -999.0;
    event.bh2mt[it]  = -999.0;
    event.bh2de[it]  = -999.0;
    event.btof[it]   = -999.0;

    event.tofmt[it]  = -999.0;
    event.tofde[it]  = -999.0;
    event.lcmt[it]   = -999.0;
    event.lcde[it]   = -999.0;
    event.stof[it]   = -999.0;
  }

  for( int it=0; it<MaxHits2; it++){
    event.tofhitpat[it] = -1;
  }
}


VEvent* ConfMan::EventAllocator()
{
  return new EventDCCheck;
}

const int NbinBcInTdc   =  500;
const double MinBcInTdc =    0.;
const double MaxBcInTdc =  500.;

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
  //***********************Chamber
  // BC1
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC1#" << std::setw(2) << i;
    title2 << "Hitpat BC1#" << std::setw(2) << i;
    title3 << "TDC First Hitpat BC1#" << std::setw(2) << i;
    title4 << "Tdc BC1#" << std::setw(2) << i;
    title5 << "Tdc First BC1#" << std::setw(2) << i;
    HB1( 100*(i+0)+0, title1.str().c_str(), MaxWireBC1+1, 0., double(MaxWireBC1+1) );
    HB1( 100*(i+0)+1, title2.str().c_str(), MaxWireBC1+1, 0., double(MaxWireBC1+1) );
    HB1( 100*(i+0)+2, title3.str().c_str(), MaxWireBC1+1, 0., double(MaxWireBC1+1) );

    HB1( 100*(i+0)+3, title4.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
    HB1( 100*(i+0)+4, title5.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
  }

  // BC2
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC2#" << std::setw(2) << i;
    title2 << "Hitpat BC2#" << std::setw(2) << i;
    title3 << "TDC First Hitpat BC2#" << std::setw(2) << i;
    title4 << "Tdc BC2#" << std::setw(2) << i;
    title5 << "Tdc First BC2#" << std::setw(2) << i;
    HB1( 100*(i+6)+0, title1.str().c_str(), MaxWireBC2+1, 0., double(MaxWireBC2+1) );
    HB1( 100*(i+6)+1, title2.str().c_str(), MaxWireBC2+1, 0., double(MaxWireBC2+1) );
    HB1( 100*(i+6)+2, title3.str().c_str(), MaxWireBC2+1, 0., double(MaxWireBC2+1) );

    HB1( 100*(i+6)+3, title4.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
    HB1( 100*(i+6)+4, title5.str().c_str(), NbinBcInTdc, MinBcInTdc, MaxBcInTdc );
  }

  // BC3
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC3#" << std::setw(2) << i;
    title2 << "Hitpat BC3#" << std::setw(2) << i;
    title3 << "TDC First Hitpat BC3#" << std::setw(2) << i;
    title4 << "Tdc BC3#" << std::setw(2) << i;
    title5 << "Tdc First BC3#" << std::setw(2) << i;
    HB1( 100*(i+12)+0, title1.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*(i+12)+1, title2.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );
    HB1( 100*(i+12)+2, title3.str().c_str(), MaxWireBC3+1, 0., double(MaxWireBC3+1) );

    HB1( 100*(i+12)+3, title4.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*(i+12)+4, title5.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
  }

  // BC4
  for( int i=1; i<=NumOfLayersBc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits BC4#" << std::setw(2) << i;
    title2 << "Hitpat BC4#" << std::setw(2) << i;
    title3 << "TDC First Hitpat BC4#" << std::setw(2) << i;
    title4 << "Tdc BC4#" << std::setw(2) << i;
    title5 << "Tdc First BC4#" << std::setw(2) << i;
    HB1( 100*(i+18)+0, title1.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+18)+1, title2.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );
    HB1( 100*(i+18)+2, title3.str().c_str(), MaxWireBC4+1, 0., double(MaxWireBC4+1) );

    HB1( 100*(i+18)+3, title4.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
    HB1( 100*(i+18)+4, title5.str().c_str(), NbinBcOutTdc, MinBcOutTdc, MaxBcOutTdc );
  }
  
  // SDC1
  for( int i=1; i<=NumOfLayersSdc-2; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC1#" << std::setw(2) << i;
    title2 << "Hitpat SDC1#" << std::setw(2) << i;
    title3 << "TDC First Hitpat SDC1#" << std::setw(2) << i;
    title4 << "Tdc SDC1#" << std::setw(2) << i;
    title5 << "Tdc First SDC1#" << std::setw(2) << i;
    HB1( 3000+100*(i+0)+0, title1.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 3000+100*(i+0)+1, title2.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );
    HB1( 3000+100*(i+0)+2, title3.str().c_str(), MaxWireSDC1+1, 0., double(MaxWireSDC1+1) );

    HB1( 3000+100*(i+0)+3, title4.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 3000+100*(i+0)+4, title5.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
  }

  // SDC2
  for( int i=1; i<=NumOfLayersSdc; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SDC2#" << std::setw(2) << i;
    title2 << "Hitpat SDC2#" << std::setw(2) << i;
    title3 << "TDC First Hitpat SDC2#" << std::setw(2) << i;
    title4 << "Tdc SDC2#" << std::setw(2) << i;
    title5 << "Tdc First SDC2#" << std::setw(2) << i;
    HB1( 3000+100*(i+4)+0, title1.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 3000+100*(i+4)+1, title2.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );
    HB1( 3000+100*(i+4)+2, title3.str().c_str(), MaxWireSDC2+1, 0., double(MaxWireSDC2+1) );

    HB1( 3000+100*(i+4)+3, title4.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
    HB1( 3000+100*(i+4)+4, title5.str().c_str(), NbinSdcInTdc, MinSdcInTdc, MaxSdcInTdc );
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

  //Tree
  HBTree("tree","tree of Counter");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  tree->Branch("trigtype",   &event.trigtype,   "trigtype/I");
  tree->Branch("trigflag",    event.trigflag,   "trigflag[10]/I");

  tree->Branch("bc1nhits", event.bc1nhits,  "bc1nhits[6]/I");
  tree->Branch("bc2nhits", event.bc2nhits,  "bc2nhits[6]/I");
  tree->Branch("bc3nhits", event.bc3nhits,  "bc3nhits[6]/I");
  tree->Branch("bc4nhits", event.bc4nhits,  "bc4nhits[6]/I");

  tree->Branch("bc1hitpat", event.bc1hitpat,  "bc1hitpat[6][100]/I");
  tree->Branch("bc2hitpat", event.bc2hitpat,  "bc2hitpat[6][100]/I");
  tree->Branch("bc3hitpat", event.bc3hitpat,  "bc3hitpat[6][100]/I");
  tree->Branch("bc4hitpat", event.bc4hitpat,  "bc4hitpat[6][100]/I");

  tree->Branch("sdc1nhits", event.sdc1nhits,  "sdc1nhits[6]/I");
  tree->Branch("sdc2nhits", event.sdc2nhits,  "sdc2nhits[6]/I");
  tree->Branch("sdc3nhits", event.sdc3nhits,  "sdc3nhits[6]/I");
  tree->Branch("sdc4nhits", event.sdc4nhits,  "sdc4nhits[6]/I");

  tree->Branch("sdc1hitpat", event.sdc1hitpat,  "sdc1hitpat[6][100]/I");
  tree->Branch("sdc2hitpat", event.sdc2hitpat,  "sdc2hitpat[6][100]/I");
  tree->Branch("sdc3hitpat", event.sdc3hitpat,  "sdc3hitpat[6][100]/I");
  tree->Branch("sdc4hitpat", event.sdc4hitpat,  "sdc4hitpat[6][100]/I");

  tree->Branch("bh1nhits",   &event.bh1nhits,   "bh1nhits/I");
  tree->Branch("bh1hitpat",   event.bh1hitpat,  "bh1hitpat[20]/I");

  tree->Branch("bh2nhits",   &event.bh2nhits,   "bh2nhits/I");
  tree->Branch("bh2hitpat",   event.bh2hitpat,  "bh2hitpat[20]/I");

  tree->Branch("tofnhits",   &event.tofnhits,   "tofnhits/I");
  tree->Branch("tofhitpat",   event.tofhitpat,  "tofhitpat[64]/I");

  tree->Branch("nhitclb",  &event.nhitclb,   "nhitclb/I");
  tree->Branch("bh1mt",     event.bh1mt,     "bh1mt[nhitclb]/D");  
  tree->Branch("bh1de",     event.bh1de,     "bh1de[nhitclb]/D");  
  tree->Branch("bh2mt",     event.bh2mt,     "bh2mt[nhitclb]/D");
  tree->Branch("bh2de",     event.bh2de,     "bh2de[nhitclb]/D");  
  tree->Branch("btof",      event.btof,      "btof[nhitclb]/D");  

  tree->Branch("nhitcls",  &event.nhitcls,   "nhitcls/I");
  tree->Branch("tofmt",     event.tofmt,     "tofmt[nhitcls]/D");  
  tree->Branch("tofde",     event.tofde,     "tofde[nhitcls]/D");  
  tree->Branch("lcmt",      event.lcmt,      "lcmt[nhitcls]/D");
  tree->Branch("lcde",      event.lcde,     " lcde[nhitcls]/D");  
  tree->Branch("stof",      event.stof,      "stof[nhitcls]/D");  

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
