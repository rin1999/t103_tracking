/*
  Hodo1Hit.hh

  202404/ K.Shirotori
*/

#ifndef Hodo1Hit_h
#define Hodo1Hit_h 1

#include "HodoRawHit.hh"
#include "ThreeVector.hh"
#include <cstddef>
#include <vector>

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
   double tot_, t_, ct_;
   double t1l_, t1h_;
   double ltdc1st_;
   double totraw_;
   
public:
   HodoRawHit * GetRawHit( void ) { return raw_; }
   bool calculate( void );
   
   bool status( void ) const { return Status_; } 
   
   double GetA( void )  const { return tot_; }
   double GetT( void )  const { return t_; }
   double GetCT( void ) const { return ct_; }
   double Get1stLtdc( void )  const { return ltdc1st_; }
   double GetTot( void )  const { return totraw_; }

   double Time( void ) const { return GetT(); }
   double DeltaE( void ) const { return GetA(); }
   
   int DetectorId( void ) const { return raw_->DetectorId(); }
   int PlaneId( void ) const { return raw_->LayerId(); }
   int SegmentId( void ) const { return raw_->SegmentId(); }
   
  virtual bool ReCalc( bool applyRecursively=false ) 
      { return calculate(); }
   
};

#endif
