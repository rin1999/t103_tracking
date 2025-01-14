/*
  Hodo1Hit.cc

  2024/04 K. Shirotori
*/

#include "Hodo1Hit.hh"

#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
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
   if(  !(raw_->GetSize_lTdc1()>0) ) return Status_;
   
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
   
   int cid=raw_->DetectorId(), layid=raw_->LayerId(), seg=raw_->SegmentId();

   //TDC
   hodoMan->GetTdcLow(cid,layid,seg,0,t1l_);
   hodoMan->GetTdcHigh(cid,layid,seg,0,t1h_);
   double tdc_cut_low = t1l_;
   double tdc_cut_High = t1h_;
   double tdc_1st = 9999.0;
   
   int size_tdc = raw_->GetSize_lTdc1();
   for( int i=0; i<size_tdc; ++i ){
      double ltdc = raw_->GetlTdc1(i);
      if( tdc_cut_low<ltdc && ltdc<tdc_cut_High ){
         if( ltdc<tdc_1st ) tdc_1st = ltdc;
      }
   }
  
   if( !hodoMan->GetTime(cid,layid,seg,0,tdc_1st,t_) ) return Status_;

   //TOT
   double tot_max = -1.0;
   int size_tot = raw_->GetSize_Tot1();
   for( int i=0; i<size_tot; ++i ){
      double tot = raw_->GetTot1(i);
      if( tot>tot_max ) tot_max = tot;
   }
   
  if( tot_max>0 ){
    if( !hodoMan->GetDe(cid,layid,seg,0,tot_max,tot_) ) return Status_;
  }
  else
    tot_=0.;

  ct_=t_;

  if( phcMan ){
    phcMan->doCorrection(cid,layid,seg,0,t_,tot_,ct_ );
  }

  return Status_=true;
}
