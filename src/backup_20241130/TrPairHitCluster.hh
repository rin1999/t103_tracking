/*
  TrPairHitCluster.hh

  2019/2  K.Shirotori
*/

#ifndef TrPairHitCluster_h
#define TrPairHitCluster_h 1

class TrLTrackHit;

class TrPairHitCluster
{
public:
  TrPairHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB=0 );
  ~TrPairHitCluster();

private:
  TrLTrackHit *hitA_, *hitB_;
  int nhits_;

public:
  int NumberOfHits( void ) const { return nhits_; }
  TrLTrackHit *GetHit( int i ) const
  {
    if(i==0)      return hitA_;
    else if(i==1) return hitB_;
    else          return 0;
  }

};

#endif
