/*
  DCPairHitCluster.cc

  2024/05 K.Shirotori
*/

#include "DCPairHitCluster.hh"

DCPairHitCluster::DCPairHitCluster( DCLTrackHit *hitA, DCLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), nhits_(0)
{
  if(hitA_) ++nhits_;
  if(hitB_) ++nhits_;
}

DCPairHitCluster::~DCPairHitCluster()
{
}
