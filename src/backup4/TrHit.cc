/*
  TrHit.cc

  2019/2  K.Shirotori
*/

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "TrHit.hh"
#include "ConfMan.hh"
#include "TrGeomMan.hh"
#include "TrTdcCalibMan.hh"
#include "TrPHCMan.hh"
#include "TrLTrackHit.hh"
#include "TrParameters.hh"

#include "DetectorID.hh"

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

TrHit::TrHit()
  : layer_(-1), fiber_(-1)
{}

TrHit::TrHit( int layer, double fiber )
  : layer_(layer), fiber_(fiber)
{}

TrHit::~TrHit()
{  
  clearRegisteredHits();
}

void TrHit::SetTdcVal( int tdc )
{
  tdc_.push_back(tdc); 
  belongTrack_.push_back(false);
  timeRange_.push_back(false);
}

void TrHit::clearRegisteredHits( void )
{
  int n=Cont_.size();
  for(int i=0; i<n; ++i)
    delete Cont_[i];
}

bool TrHit::CalcTrObservables( void )
{
  static const std::string funcname="[TrHit::CalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  TrGeomMan *geomMan=confMan->GetTrGeomManager();
  if(!geomMan) return false;
  TrTdcCalibMan *calibMan=confMan->GetTrTdcCalibManager();
  if(!calibMan) return false;
  TrPHCMan *phcMan=confMan->GetTrPHCManager();
  if(!phcMan) return false;

  fpos_ =geomMan->calcWirePosition(layer_,fiber_);
  mfpos_=geomMan->calcWirePosition(layer_,mfiber_);
  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhittdc  = tdc_.size();
  for (int i=0; i<nhittdc; i++) {
    double ctime;
    if(!calibMan->GetTime( layer_, fiber_, tdc_[i], ctime ))
      return false;
    
    time_.push_back(ctime);
    dl_.push_back(0.0);

    //bool status=phcMan->Correction( layer_, wire_, ctime, dtime, dlength );
    //if (status == false) Status = status;
    
    if( ctime>MinTimeTr[layer_-PlOffsSFT] && 
	ctime<MaxTimeTr[layer_-PlOffsSFT] ){
      timeRange_[i]=true;
    }
    else Status = false;
    
    if( Status ){
      std::cout<< layer_ << "-" 
      	       << fiber_ << " : " 
      	       << tdc_[i] << " -> " << ctime << std::endl;    
    }
  }

  return Status;
}

bool TrHit::ReCalcTrObservables( void )
{
  static const std::string funcname="[TrHit::ReCalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  TrGeomMan *geomMan=confMan->GetTrGeomManager();
  if(!geomMan) return false;

  mfpos_=geomMan->calcWirePosition(layer_,mfiber_);
  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhittime  = mtime_.size();
  for (int i=0; i<nhittime; i++) {
    double ctime = mtime_[i];

    dl_.push_back(0.0);

    if( ctime>MinTimeTr[layer_-PlOffsSFT] && 
    	ctime<MaxTimeTr[layer_-PlOffsSFT] ){
      timeRange_[i]=true;
    }
    else Status = false;
    
    if( Status ){
      // std::cout<< layer_ << "-" 
      // 	       << mfiber_ << " : " 
      // 	       << mtime_[i] << std::endl;    
    }
  }

  return Status;
}

