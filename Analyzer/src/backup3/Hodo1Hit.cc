/*
  Hodo1Hit.cc
*/

#include "Hodo1Hit.hh"

#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
//#include "SksObjectId.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdexcept>
#include <sstream>


Hodo1Hit::Hodo1Hit( HodoRawHit *rhit )
  : raw_(rhit), Status_(false)
{}

Hodo1Hit::~Hodo1Hit()
{}

bool Hodo1Hit::calculate( void )
{
  static const std::string funcname = "[Hodo1Hit::calculate]";

  Status_=false;
  if( raw_->GetNumOfTdcHits()!=1 ) return Status_;

  int tdc=raw_->GetTdc1();
  int adc=raw_->GetAdc1();
  int UorD=0;
  if( tdc<0 ){
    tdc=raw_->GetTdc2();  
    adc=raw_->GetAdc2();  
    UorD=1;
  }

  if( tdc<0 ) return Status_;

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << funcname << ": cannot get confManager" << std::endl; 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  HodoPHCMan   *phcMan  = confMan->GetHodoPHCManager();
  if(!hodoMan){
    std::cerr << funcname << ": cannot get HodoParamManager" << std::endl; 
    return Status_;
  }

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(),
    seg=raw_->SegmentId();

  if( !hodoMan->GetTime(cid,plid-1,seg,UorD,tdc,t_) ) return Status_;

  if( adc>=0 ){
    if( !hodoMan->GetDe(cid,plid-1,seg,UorD,adc,a_) ) return Status_;
  }
  else
    a_=0.;

  ct_=t_;
  if( phcMan ){
    phcMan->doCorrection(cid,plid-1,seg,UorD,t_,a_,ct_ );
  }

  return Status_=true;
}
