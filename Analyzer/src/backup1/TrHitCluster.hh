/*
  TrHitCluster.hh

  2012/5 K.Shirotori
*/

#ifndef TrHitCluster_h
#define TrHitCluster_h 1

class TrLTrackHit;

class TrHitCluster
{
public:
  TrHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB=0 );
  ~TrHitCluster();

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
