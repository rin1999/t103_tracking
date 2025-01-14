/*
  Hodo2Hit.cc

  2024/04 K. Shirotori
*/

#include "Hodo2Hit.hh"

#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

Hodo2Hit::Hodo2Hit( HodoRawHit *rhit )
   : raw_(rhit), Status_(false)
{
}

Hodo2Hit::~Hodo2Hit()
{
}

double Hodo2Hit::DeltaE( void ) const
{
   return sqrt(fabs(tot1_*tot2_));
}

bool Hodo2Hit::calculate( void )
{
   static const std::string funcname = "[Hodo2Hit::calculate]";

   Status_=false;
   if( !(raw_->GetSize_lTdc1()>0 && raw_->GetSize_lTdc2()>0) ) return Status_;
    
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

   //std::cout<< "***********" << std::endl;
   int cid=raw_->DetectorId(), layid=raw_->LayerId(), seg=raw_->SegmentId();
   //std::cout<< cid << "-" << layid << "-" << seg <<std::endl;
   
   //TDC
   if( !hodoMan->GetTdcLow(cid,layid,seg,0,t1l_) ||
       !hodoMan->GetTdcHigh(cid,layid,seg,0,t1h_) ||
       !hodoMan->GetTdcLow(cid,layid,seg,1,t2l_) ||
       !hodoMan->GetTdcHigh(cid,layid,seg,1,t2h_) )  return Status_;

   //std::cout<< t1l_ << ":" << t1h_ << "-" << t2l_ << ":" << t1h_ << std::endl;
 
   double tdc1_cut_low  = t1l_; 
   double tdc1_cut_high = t1h_;
   double tdc2_cut_low  = t2l_; 
   double tdc2_cut_high = t2h_;

   double tdc1_1st = 9999.0;
   double tdc2_1st = 9999.0;
   
   int size_tdc1 = raw_->GetSize_lTdc1();
   for( int i1=0; i1<size_tdc1; ++i1 ){
      double ltdc1 = raw_->GetlTdc1(i1);
      if( tdc1_cut_low<ltdc1 && ltdc1<tdc1_cut_high ){
         if( ltdc1<tdc1_1st ) tdc1_1st = ltdc1;
      }
      else tdc1_1st = -1.0;
   }
   int size_tdc2=raw_->GetSize_lTdc2();
   for( int i2=0; i2<size_tdc2; ++i2 ){
      double ltdc2 = raw_->GetlTdc2(i2);
      if( tdc2_cut_low<ltdc2 && ltdc2<tdc2_cut_high ){
         if( ltdc2<tdc2_1st ) tdc2_1st = ltdc2;
      }
      else tdc2_1st = -1.0;
   }

   // std::cout<< "******" << std::endl;
   // std::cout<< tdc1_1st << " : " << tdc2_1st << std::endl;

   if( !(tdc1_1st>0 && tdc2_1st>0) ) return Status_;
   
   if( !hodoMan->GetTime(cid,layid,seg,0,tdc1_1st,t1_) ||
       !hodoMan->GetTime(cid,layid,seg,1,tdc2_1st,t2_) ) return Status_;

   // std::cout<< t1_ <<" : " << t2_ << std::endl;

   //TOT
   double tot1_max = -1.0;
   double tot2_max = -1.0;
   
   int size_tot1 = raw_->GetSize_Tot1();
   for( int i=0; i<size_tot1; ++i ){
      double tot1 = raw_->GetTot1(i);
      if( tot1>tot1_max ) tot1_max = tot1;
   }
   int size_tot2 = raw_->GetSize_Tot2();
   for( int i=0; i<size_tot2; ++i ){
      double tot2 = raw_->GetTot2(i);
      if( tot2>tot2_max ) tot2_max = tot2;
   }

   if( tot1_max>0 ){
      if( !hodoMan->GetDe(cid,layid,seg,0,tot1_max,tot1_) ) return Status_;
   }
   else
      tot1_=0.;

   if( tot2_max>0 ){
      if( !hodoMan->GetDe(cid,layid,seg,1,tot2_max,tot2_) ) return Status_;
   }
   else
      tot2_=0.;
   
   ct1_=t1_; ct2_=t2_;
   
  // if( phcMan ){
  //   phcMan->doCorrection(cid,layid,seg,0,t1_,tot1_,ct1_ );
  //   phcMan->doCorrection(cid,layid,seg,1,t2_,tot2_,ct2_ );
  // }

  return Status_=true;
} 
