/*
  TrHitCluster.cc

  2012/5 K.Shirotori
*/

#include "TrHitCluster.hh"

TrHitCluster::TrHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), nhits_(0)
{
  if(hitA_) ++nhits_;
  if(hitB_) ++nhits_;
}

TrHitCluster::~TrHitCluster()
{
}
