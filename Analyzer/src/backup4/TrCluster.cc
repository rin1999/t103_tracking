/*
  TrCluster.cc

  2018/2  K.Shirotori
*/

#include "TrCluster.hh"
#include "TrHit.hh"
#include "TrAnalyzer.hh"

#include <cstring>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

TrCluster::TrCluster( TrHit *hitA, TrHit *hitB, TrHit *hitC )
  : hitA_(hitA), hitB_(hitB), hitC_(hitC), csize_(0), gfastatus_(true)
{
  if(hitA) ++csize_;
  if(hitB) ++csize_;
  if(hitC) ++csize_;

  calculate();
}

TrHit * TrCluster::GetHit( int i ) const
{
  if(i==0) return hitA_;
  else if(i==1) return hitB_;
  else if(i==2) return hitC_;
  else return 0;
}

void TrCluster::calculate( void )
{
  double mf=0., mt=0.;
  if( hitA_ ){
    mf += hitA_->GetFiber();
    mt += hitA_->GetTime(0);
  }
  if( hitB_ ){
    mf += hitB_->GetFiber();
    mt += hitB_->GetTime(0);
  }
  if( hitC_ ){
    mf += hitC_->GetFiber();
    mt += hitC_->GetTime(0);
  }
  mf /= double(csize_);
  mt /= double(csize_);

  MeanFiber_=mf;
  MeanTime_=mt; 

  if(hitA_){
    hitA_->SetMeanFiber(mf);
    hitA_->SetMeanTime(mt);
  }
  if(hitB_){
    hitB_->SetMeanFiber(mf);
    hitB_->SetMeanTime(mt);
  }
  if(hitC_){
    hitC_->SetMeanFiber(mf);
    hitC_->SetMeanTime(mt);
  }
} 

bool TrCluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    if(hitA_) hitA_->ReCalcTr(applyRecursively);
    if(hitB_) hitB_->ReCalcTr(applyRecursively);
    if(hitC_) hitC_->ReCalcTr(applyRecursively);
  }
  calculate();

  return true;
}
