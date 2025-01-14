/*
  TrPairHitCluster.cc

  2019/2  K.Shirotori
*/

#include "TrPairHitCluster.hh"

TrPairHitCluster::TrPairHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), nhits_(0)
{
  if(hitA_) ++nhits_;
  if(hitB_) ++nhits_;
}

TrPairHitCluster::~TrPairHitCluster()
{
}
