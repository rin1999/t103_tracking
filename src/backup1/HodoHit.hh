/*
  HodoHit.hh
*/

#ifndef HodoHit_h
#define HodoHit_h 1

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;

class RawData;

class HodoHit
{
public:
  explicit HodoHit( HodoRawHit *rhit );
  virtual ~HodoHit();

private:
  HodoHit( const HodoHit & );
  HodoHit & operator = ( const HodoHit & );

protected:
  HodoRawHit *raw_;
  bool Status_;
  DoubleVec PosX_, PosY_;
  DoubleVec Time_, Edep_;
  DoubleVec Mom_, Beta_;

public:
  HodoRawHit * GetRawHit( void ) { return raw_; }
  bool calculate( void );

  bool status( void ) const { return Status_; } 

  int DetectorId( void ) const { return raw_->DetectorId(); }
  int PlaneId( void ) const { return raw_->LayerId(); }
  int SegmentId( void ) const { return raw_->SegmentId(); }

  double GetTime( int nh ) const { return raw_->GetTime(nh); }
  double GetEdep( int nh ) const { return raw_->GetEdep(nh); }
  int    GetSize( void ) const { return raw_->GetSize(); }

  virtual bool ReCalc( bool applyRecursively=false ) 
  { return calculate(); }

};

#endif
