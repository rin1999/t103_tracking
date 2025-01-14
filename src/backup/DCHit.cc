/*
  DCHit.cc

  2024/05  K.Shirotori
*/

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "DCHit.hh"
#include "ConfMan.hh"
#include "GeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "DCLTrackHit.hh"
#include "DCParameters.hh"

#include "DetectorInfo.hh"

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

DCHit::DCHit()
  : layer_(-1), wire_(-1)
{}

DCHit::DCHit( int layer, double wire )
  : layer_(layer), wire_(wire)
{}

DCHit::~DCHit()
{  
  clearRegisteredHits();
}

void DCHit::SetTdcVal( double tdc )
{
  tdc_.push_back(tdc); 
  belongTrack_.push_back(false);
  dlRange_.push_back(false);
}

void DCHit::clearRegisteredHits( void )
{
  int n=Cont_.size();
  for(int i=0; i<n; ++i)
    delete Cont_[i];
}

bool DCHit::CalcDCObservables( void )
{
  static const std::string funcname="[DCHit::CalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  GeomMan *geomMan=confMan->GetGeomManager();
  if(!geomMan) return false;
  DCTdcCalibMan *calibMan=confMan->GetDCTdcCalibManager();
  if(!calibMan) return false;
  DCDriftParamMan *driftMan=confMan->GetDCDriftParamManager();
  if(!driftMan) return false;

  wpos_=geomMan->calcWirePosition(layer_,wire_);
  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhitdc  = tdc_.size();
  for (int i=0; i<nhitdc; i++) {
    double ctime;
    //std::cout<< tdc_[i] << std::endl;
    if(!calibMan->GetTime( layer_, wire_, tdc_[i], ctime ))
      return false;

    double dtime, dlength;
    bool status=driftMan->calcDrift( layer_, wire_, ctime, dtime, dlength );

    // std::cout<< layer_ << " : " << wire_ << " : "
    //          <<tdc_[i] << " : " << dtime << " : " << dlength
    //          << std::endl;
    
    if (status == false) Status = status;
    
    dt_.push_back(dtime);
    dl_.push_back(dlength);
    
    if(layer_>=101 && layer_<=108 ){
       if( dl_[i]>MinDLBDC[layer_-PlOffsBDC] && dl_[i]<MaxDLBDC[layer_-PlOffsBDC] )
          dlRange_[i]=true;
    }
    if(layer_>=201 && layer_<=208 ){
       if( dl_[i]>MinDLKLDC[layer_-PlOffsKLDC] && dl_[i]<MaxDLKLDC[layer_-PlOffsKLDC] )
          dlRange_[i]=true;
    }
  }
  
  return Status;
}

