/*
  UserDCHitCheck.cc
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "SksLib.hh"
#include "DetectorID.hh"
#include "HistHelper.hh"

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventDCHitCheck : public VEvent
{

private:
  RawData *rawData;
  DCAnalyzer *DCAna;

public:
  EventDCHitCheck();
  ~EventDCHitCheck();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventDCHitCheck::EventDCHitCheck()
  : VEvent(),
    rawData(0),
    DCAna(new DCAnalyzer)
{
}

EventDCHitCheck::~EventDCHitCheck()
{
  delete DCAna;
  if (rawData) delete rawData;
}

bool EventDCHitCheck::ProcessingBegin()
{
  //  ConfMan::InitializeHistograms(); 
 return true;
}

bool EventDCHitCheck::ProcessingNormal( void )
{
  static const std::string funcname = 
    "[EventDCHitCheck::ProcessingNormal]";

  rawData = new RawData;
  rawData->DecodeHits();

  DCAna->DecodeRawHits( rawData );

  // BCIn
  for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
    int pid=layer+PlMinBcIn-1;
    const DCHitContainer & cont = DCAna->GetBcInHC(layer);
    int nh=cont.size();
#if 0
    std::cout << funcname << ": BcIn " << layer 
	      << " " << nh << std::endl;
#endif
    HF1( 100*pid, double(nh) );
    for( int i=0; i<nh; ++i ){
      DCHit *hit=cont[i];
      double wire=hit->GetWire();

      HF1(100*pid+1, wire-0.5 );
      int nhdt = hit->GetDriftTimeSize();
      for( int k=0; k<nhdt; k++ ){
	double dt=hit->GetDriftTime(k);
	HF1(100*pid+2, dt );
      }
      int nhdl = hit->GetDriftLengthSize();
      for( int k=0; k<nhdl; k++ ){
	double dl=hit->GetDriftLength(k);
	HF1(100*pid+3, dl );
      }
    }
  }
  // BCOut
  for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
    int pid=layer+PlMinBcOut-1;
    const DCHitContainer & cont = DCAna->GetBcOutHC(layer);
    int nh=cont.size();
#if 0
    std::cout << funcname << ": BcOut " << layer 
	      << " " << nh << std::endl;
#endif
    HF1( 100*pid, double(nh) );
    for( int i=0; i<nh; ++i ){
      DCHit *hit=cont[i];
      double wire=hit->GetWire();

      HF1(100*pid+1, wire-0.5 );
      int nhdt = hit->GetDriftTimeSize();
      for( int k=0; k<nhdt; k++ ){
	double dt=hit->GetDriftTime(k);
	HF1(100*pid+2, dt );
      }
      int nhdl = hit->GetDriftLengthSize();
      for( int k=0; k<nhdl; k++ ){
	double dl=hit->GetDriftLength(k);
	HF1(100*pid+3, dl );
      }
    }
  }

  // SDCIn
  for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
    int pid=layer+PlMinSdcIn-1+30;
    const DCHitContainer & cont = DCAna->GetSdcInHC(layer);
    int nh=cont.size();
#if 0
    std::cout << funcname << ": SdcIn " << layer 
	      << " " << nh << std::endl;
#endif
    HF1( 100*pid, double(nh) );
    for( int i=0; i<nh; ++i ){
      DCHit *hit=cont[i];
      double wire=hit->GetWire();

      HF1(100*pid+1, wire-0.5 );
      int nhdt = hit->GetDriftTimeSize();
      for( int k=0; k<nhdt; k++ ){
	double dt=hit->GetDriftTime(k);
	HF1(100*pid+2, dt );
      }
      int nhdl = hit->GetDriftLengthSize();
      for( int k=0; k<nhdl; k++ ){
	double dl=hit->GetDriftLength(k);
	HF1(100*pid+3, dl );
      }
    }
  }
    
  // SDCOut
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    int pid=layer+PlMinSdcOut-1+20;
    const DCHitContainer & cont = DCAna->GetSdcOutHC(layer);
    int nh=cont.size();
#if 0
    std::cout << funcname << ": SdcOut " << layer 
	      << " " << nh << std::endl;
#endif
    HF1( 100*pid, double(nh) );
    for( int i=0; i<nh; ++i ){
      DCHit *hit=cont[i];
      double wire=hit->GetWire();

      HF1(100*pid+1, wire-0.5 );
      int nhdt = hit->GetDriftTimeSize();
      for( int k=0; k<nhdt; k++ ){
	double dt=hit->GetDriftTime(k);
	HF1(100*pid+2, dt );
      }
      int nhdl = hit->GetDriftLengthSize();
      for( int k=0; k<nhdl; k++ ){
	double dl=hit->GetDriftLength(k);
	HF1(100*pid+3, dl );
      }
    }
  }

  return true;

}

bool EventDCHitCheck::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent* ConfMan::EventAllocator()
{
  return new EventDCHitCheck;
}


const int NBin1DTBc = 500;
const double MinDTBc = -100.;
const double MaxDTBc =  400.;
const int NBin1DLBc = 500;
const double MinDLBc = -1.0;
const double MaxDLBc =  4.0;

const int NBin1DTSdcIn = 500;
const double MinDTSdcIn = -100.;
const double MaxDTSdcIn =  400.;
const int NBin1DLSdcIn = 500;
const double MinDLSdcIn = -1.0;
const double MaxDLSdcIn =  4.0;


const int NBin1DTSdcOut = 700;
const double MinDTSdcOut = -100.;
const double MaxDTSdcOut =  600.;
const int NBin1DLSdcOut =  800;
const double MinDLSdcOut = -5.0;
const double MaxDLSdcOut = 35.0;

bool ConfMan::InitializeHistograms( void )
{
  // BC
  for( int i=1; i<=24; ++i ){
    std::ostringstream title0, title1, title2, title3;
    title0 << "#Hits Bc" << std::setw(2) << i;
    HB1( 100*i, title0.str().c_str(), 49, 0., 49. );
    title1 << "HitPat Bc" << std::setw(2) << i;
    HB1( 100*i+1, title1.str().c_str(), 48, 0., 48. );
    title2 << "DriftTime Bc" << std::setw(2) << i;
    HB1( 100*i+2, title2.str().c_str(), NBin1DTBc, MinDTBc, MaxDTBc );
    title3 << "DriftLength Bc" << std::setw(2) << i;
    HB1( 100*i+3, title3.str().c_str(), NBin1DLBc, MinDLBc, MaxDLBc );
  }

  // SDC1
  for( int i=1; i<=5; ++i ){
    std::ostringstream title0, title1, title2, title3;
    title0 << "#Hits Sdc" << std::setw(2) << i;
    HB1( 100*i+3000, title0.str().c_str(), 129, 0., 129. );
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+3001, title1.str().c_str(), 128, 0., 128. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+3002, title2.str().c_str(), 
	 NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+3003, title3.str().c_str(), 
	 NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
  }

  // SDC2
  for( int i=6; i<=11; ++i ){
    std::ostringstream title0, title1, title2, title3;
    title0 << "#Hits Sdc" << std::setw(2) << i;
    HB1( 100*i+3000, title0.str().c_str(), 81, 0., 81. );
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+3001, title1.str().c_str(), 80, 0., 80. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+3002, title2.str().c_str(), 
	 NBin1DTSdcIn, MinDTSdcIn, MaxDTSdcIn );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+3003, title3.str().c_str(), 
	 NBin1DLSdcIn, MinDLSdcIn, MaxDLSdcIn );
  }

  // SDC34
  for( int i=31; i<=46; ++i ){
    std::ostringstream title0, title1, title2, title3;
    title0 << "#Hits Sdc" << std::setw(2) << i;
    HB1( 100*i+2000, title0.str().c_str(), 25, 0., 25. );
    title1 << "HitPat Sdc" << std::setw(2) << i;
    HB1( 100*i+2001, title1.str().c_str(), 24, 0., 24. );
    title2 << "DriftTime Sdc" << std::setw(2) << i;
    HB1( 100*i+2002, title2.str().c_str(), 
	 NBin1DTSdcOut, MinDTSdcOut, MaxDTSdcOut );
    title3 << "DriftLength Sdc" << std::setw(2) << i;
    HB1( 100*i+2003, title3.str().c_str(), 
	 NBin1DLSdcOut, MinDLSdcOut, MaxDLSdcOut );
  }
    
  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
//   if( CMapFileName_!="" )
//     CMapManager_ = new CMapMan(CMapFileName_);
//   if(CMapManager_) CMapManager_->Initialize();

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

//   if( ScalerDefinitionFileName_!="" )
//     ScalerAnalyzer_ = new ScalerAna(ScalerDefinitionFileName_);
//   else
//     ScalerAnalyzer_ = new ScalerAna();
//   if(ScalerAnalyzer_) ScalerAnalyzer_->Initialize(); 

  return true;
}
