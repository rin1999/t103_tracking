/*
  DCHit.cc

  2018/12  K.Shirotori
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
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "DCLTrackHit.hh"
#include "DCParameters.hh"

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

void DCHit::SetTdcVal( int tdc )
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
  DCGeomMan *geomMan=confMan->GetDCGeomManager();
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
    
    if (status == false) Status = status;
    
    dt_.push_back(dtime);
    dl_.push_back(dlength);
    
    if(layer_>=200){
      if( dl_[i]>MinDLDC[layer_-200] && dl_[i]<MaxDLDC[layer_-200] )
    dlRange_[i]=true;
    }
  }

  return Status;
}

