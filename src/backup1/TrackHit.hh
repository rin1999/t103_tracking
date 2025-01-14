/*
  TrackHit.hh

  2012/5  K.Shirotori
*/


#ifndef TrackHit_h
#define TrackHit_h 1

#include "TrLTrackHit.hh" 
#include "ThreeVector.hh"

#include <cstddef>

class TrLocalTrack;

class TrackHit
{
public:
  explicit TrackHit( TrLTrackHit *hit );
  ~TrackHit();
private:
  TrackHit( const TrackHit & );
  TrackHit & operator= ( const TrackHit & );

public:
  TrLTrackHit *GetHit( void ) { return trhitp_; }
  void SetCalGPos( const ThreeVector &pos ) { calGPos_=pos; }
  void SetCalLPos( double pos ) { calLPos_=pos; }
  
  int GetLayer( void ) const { return trhitp_->GetLayer(); }
  double GetLocalHitPos( void ) const 
  { return trhitp_->GetLocalHitPos(); }
  const ThreeVector & GetCalGPos( void ) const 
  { return calGPos_; }
  double GetCalLPos( void ) const { return calLPos_; }
  double GetResidual( void ) const 
  { return trhitp_->GetLocalHitPos()-calLPos_; }

  bool ReCalc( bool applyRecursively=false )
  {
    if( applyRecursively )
      return trhitp_->ReCalc(applyRecursively);
    else
      return true;
  }

private:
  TrLTrackHit *trhitp_;
  ThreeVector calGPos_;
  double calLPos_;
};

#endif

