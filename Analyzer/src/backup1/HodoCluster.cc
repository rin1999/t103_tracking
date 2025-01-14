/*
  HodoCluster.cc
*/

#include "HodoCluster.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"

#include <cstring>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

HodoCluster::HodoCluster( HodoHit *hitA, HodoHit *hitB, HodoHit *hitC )
  : hitA_(hitA), hitB_(hitB), hitC_(hitC), csize_(0), gfastatus_(true)
{
  if(hitA) ++csize_;
  if(hitB) ++csize_;
  if(hitC) ++csize_;

  calculate();
}

HodoHit * HodoCluster::GetHit( int i ) const
{
  if(i==0) return hitA_;
  else if(i==1) return hitB_;
  else if(i==2) return hitC_;
  else return 0;
}

void HodoCluster::calculate( void )
{
  double ms=0., mt=0., de=0., dt=0.;
  if( hitA_ ){
    ms += hitA_->SegmentId();
    mt += hitA_->GetTime(0);
    de += hitA_->GetEdep(0);
    dt += (hitA_->GetTime(0));
  }
  if( hitB_ ){
    ms += hitB_->SegmentId();
    mt += hitB_->GetTime(0);
    de += hitB_->GetEdep(0);
    dt += (hitB_->GetTime(0));
  }
  if( hitC_ ){
    ms += hitC_->SegmentId();
    mt += hitC_->GetTime(0);
    de += hitC_->GetEdep(0);
    dt += (hitC_->GetTime(0));
  }
  ms /= double(csize_);
  mt /= double(csize_);
  dt /= double(csize_);

  MeanTime_=mt; dE_=de; MeanSeg_=ms, Tdiff_=dt;
} 

bool HodoCluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    if(hitA_) hitA_->ReCalc(applyRecursively);
    if(hitB_) hitB_->ReCalc(applyRecursively);
    if(hitC_) hitC_->ReCalc(applyRecursively);
  }
  calculate();

  return true;
}
