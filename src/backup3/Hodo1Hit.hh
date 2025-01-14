/*
  Hodo1Hit.hh
*/

#ifndef Hodo1Hit_h
#define Hodo1Hit_h 1

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

class RawData;

class Hodo1Hit
{
public:
  explicit Hodo1Hit( HodoRawHit *rhit );
  virtual ~Hodo1Hit();

private:
  Hodo1Hit( const Hodo1Hit & );
  Hodo1Hit & operator = ( const Hodo1Hit & );

protected:
  HodoRawHit *raw_;
  bool Status_;
  double a_, t_, ct_;

public:
  HodoRawHit * GetRawHit( void ) { return raw_; }
  bool calculate( void );

  bool status( void ) const { return Status_; } 

  double GetA( void )  const { return a_; }
  double GetT( void )  const { return t_; }
  double GetCT( void ) const { return ct_; }

  double Time( void ) const { return GetT(); }
  double DeltaE( void ) const { return GetA(); }

  int DetectorId( void ) const { return raw_->DetectorId(); }
  int PlaneId( void ) const { return raw_->PlaneId(); }
  int SegmentId( void ) const { return raw_->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false ) 
  { return calculate(); }

};

#endif
