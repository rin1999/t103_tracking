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

const double TdcLow  =  700.;
const double TdcHigh = 1000.;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventDCMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventDCMonitor();
  ~EventDCMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
};

EventDCMonitor::EventDCMonitor()
  : VEvent(),
    rawData(0)
{
}

EventDCMonitor::~EventDCMonitor()
{
  if (rawData) delete rawData;
}

bool EventDCMonitor::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventDCMonitor::ProcessingNormal()
{
  rawData = new RawData;
  rawData->DecodeHits();

  //**************************************************************************
  //******************RawData

  //BC1&BC2
  {
    for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
      if( layer<=NumOfLayersBc ){
	const DCRHitContainer &contIn =rawData->GetBcInRawHC(layer);
	int nhIn=contIn.size();
	HF1( 100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 100*layer+1, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( 0< tdc && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 100*layer+2, wire-0.5 );
	  HF1( 100*layer+4, tdc1st );
	}
      }
      else{
	const DCRHitContainer &contIn =rawData->GetBcInRawHC(layer);
	int nhIn=contIn.size();
	HF1( 100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 100*layer+1, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( 0< tdc && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 100*layer+2, wire-0.5 );
	  HF1( 100*layer+4, tdc1st );
	}
      }
    }
  }  

  //BC3&BC4
  int bcIn_planeNum=12;
  {
    for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
      if( layer<=NumOfLayersBc ){
	const DCRHitContainer &contOut =rawData->GetBcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 100*(layer+bcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+bcIn_planeNum)+1, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*(layer+bcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 100*(layer+bcIn_planeNum)+2, wire-0.5 );
	  HF1( 100*(layer+bcIn_planeNum)+4, tdc1st );
	}
      }
      else{
	const DCRHitContainer &contOut =rawData->GetBcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 100*(layer+bcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 100*(layer+bcIn_planeNum)+1, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 100*(layer+bcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 100*(layer+bcIn_planeNum)+2, wire-0.5 );
	  HF1( 100*(layer+bcIn_planeNum)+4, tdc1st );
	}
      }
    }
  }  

  //SDC1&SDC2
  {
    for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
      if( layer<=NumOfLayersSdc-2 ){
	const DCRHitContainer &contIn =rawData->GetSdcInRawHC(layer);
	int nhIn=contIn.size();
	HF1( 3000+100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*layer+1, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*layer+2, wire-0.5 );
	  HF1( 3000+100*layer+4, tdc1st );
	}
      }
      else{
	const DCRHitContainer &contIn =rawData->GetSdcInRawHC(layer);
	int nhIn=contIn.size();
	HF1( 3000+100*layer+0, nhIn );
	for( int i=0; i<nhIn; ++i ){
	  DCRawHit *hit=contIn[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*layer+1, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( (TdcLow< tdc && tdc < TdcHigh) && tdc > tdc1st ) tdc1st=tdc;
	    HF1( 3000+100*layer+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*layer+2, wire-0.5 );
	  HF1( 3000+100*layer+4, tdc1st );
	}
      }
    }
  }

  //SDC3&SDC4
  {
    int sdcIn_planeNum=10;
    for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
      if( layer<=NumOfLayersSdc ){
	const DCRHitContainer &contOut =rawData->GetSdcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 3000+100*(layer+sdcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+sdcIn_planeNum)+1, wire-0.5 );
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for( int k=0; k<nhtdc; k++ ){
	    int tdc = hit->GetTdc(k);
	    if( tdc ) tdc1st=tdc;
	    HF1( 3000+100*(layer+sdcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*(layer+sdcIn_planeNum)+2, wire-0.5 );
	  HF1( 3000+100*(layer+sdcIn_planeNum)+4, tdc1st );
	}
      }
      else{
	const DCRHitContainer &contOut =rawData->GetSdcOutRawHC(layer);
	int nhOut=contOut.size();
	HF1( 3000+100*(layer+sdcIn_planeNum)+0, nhOut );
	for( int i=0; i<nhOut; ++i ){
	  DCRawHit *hit=contOut[i];
	  int wire=hit->WireId();

	  HF1( 3000+100*(layer+sdcIn_planeNum)+1, wire-0.5 );	  
	  int nhtdc = hit->GetTdcSize();
	  int tdc1st = -1;
	  for(int k=0; k<nhtdc; k++){
	    int tdc = hit->GetTdc(k);
	    if( tdc ) tdc1st=tdc;
	    HF1( 3000+100*(layer+sdcIn_planeNum)+3, tdc );
	  }
	  if( tdc1st > 0 ) HF1( 3000+100*(layer+sdcIn_planeNum)+2, wire-0.5 );
	  HF1( 3000+100*(layer+sdcIn_planeNum)+4, tdc1st );
	}
      }
    }
  }

  return true;
}

bool EventDCMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventDCMonitor;
}

const int NbinBcInTdc   =  50;
const double MinBcInTdc =    0.;
const double MaxBcInTdc =  50.;

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
