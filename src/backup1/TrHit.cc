/*
  TrHit.cc

  2012/5  K.Shirotori
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
// #include "TrTdcCalibMan.hh"
// #include "TrDriftParamMan.hh"
#include "TrLTrackHit.hh"
#include "TrParameters.hh"

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

TrHit::TrHit()
  : layer_(-1), wire_(-1)
{
}

TrHit::TrHit( int layer, double wire )
  : layer_(layer), wire_(wire)
{
}

TrHit::~TrHit()
{
  clearRegisteredHits();
}

void TrHit::SetPos( double pos )
{
  pos_.push_back(pos); 
  belongTrack_.push_back(false);
  dlRange_.push_back(true);
}

void TrHit::clearRegisteredHits( void )
{
  int n=Cont_.size();
  for(int i=0; i<n; ++i)
    delete Cont_[i];
}

bool TrHit::CalcObservables( void )
{
  static const std::string funcname="[TrHit::CalcObservables]";
  
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;

  TrGeomMan *geomMan=confMan->GetTrGeomManager();
  if(!geomMan) return false;
//   TrTdcCalibMan *calibMan=confMan->GetTrTdcCalibManager();
//   if(!calibMan) return false;
//   TrDriftParamMan *driftMan=confMan->GetTrDriftParamManager();
//   if(!driftMan) return false;

  wpos_=geomMan->calcWirePosition(layer_,wire_);
  angle_=geomMan->GetTiltAngle(layer_);
  
  bool Status = true;
  int nhpos  = pos_.size();
  for (int i=0; i<nhpos; i++) {
    //     double ctime;
    //     if(!calibMan->GetTime( layer_, wire_, tdc_[i], ctime ))
    //       return false;
    
//     double dtime, dlength;
//     bool status=driftMan->calcDrift( layer_, wire_, ctime, dtime, dlength );
    
//     if (status == false) Status = status;
    
//     dt_.push_back(dtime);
    dl_.push_back(pos_[i]);
    
//     if(layer_>=100){
//       if( dl_[i]>MinDLBc[layer_-100] && dl_[i]<MaxDLBc[layer_-100] )
// 	dlRange_[i]=true;
//     }
//     else{
//       if( dl_[i]>MinDLSdc[layer_] && dl_[i]<MaxDLSdc[layer_] )
// 	dlRange_[i]=true;
//     }
//   }
  }

  return Status;
}
