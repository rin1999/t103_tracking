/*
  Hodo2Hit.hh

  2024/04 K. Shirotori
*/

#ifndef HODO2_HIT_H 
#define HODO2_HIT_H
 
#include "HodoRawHit.hh"
#include "ThreeVector.hh"
#include <cstddef>
#include <vector>

class RawData;

class Hodo2Hit
{
public:
   explicit Hodo2Hit( HodoRawHit *rhit );
   virtual ~Hodo2Hit();
   
private:
   Hodo2Hit( const Hodo2Hit & );
   Hodo2Hit & operator = ( const Hodo2Hit & );
   
protected:
   HodoRawHit *raw_;
   bool Status_;
   double t1_, t2_, ct1_, ct2_, tot1_, tot2_;
   double t1l_, t2l_, t1h_, t2h_;
   double ltdcraw_1st_1_, ltdcraw_1st_2_;
   double totraw_1_, totraw_2_;
   
public:
   HodoRawHit * GetRawHit( void ) { return raw_; }
   bool calculate( void );
   
   bool status( void ) const { return Status_; } 

   double GetTUp( void )    const { return t1_; }
   double GetTLeft( void )  const { return t1_; }
   double GetTDown( void )  const { return t2_; }
   double GetTRight( void ) const { return t2_; }

   double Get1stLtdcUp( void )    const { return ltdcraw_1st_1_; }
   double Get1stLtdcLeft( void )  const { return ltdcraw_1st_1_; }
   double Get1stLtdcDown( void )  const { return ltdcraw_1st_2_; }
   double Get1stLtdcRight( void ) const { return ltdcraw_1st_2_; }

   double GetTotUp( void )    const { return totraw_1_; }
   double GetTotLeft( void )  const { return totraw_1_; }
   double GetTotDown( void )  const { return totraw_2_; }
   double GetTotRight( void ) const { return totraw_2_; }
   
   double GetCTUp( void )    const { return ct1_; }
   double GetCTLeft( void )  const { return ct1_; }
   double GetCTDown( void )  const { return ct2_; }
   double GetCTRight( void ) const { return ct2_; }
   
   double GetAUp( void )    const { return tot1_; }
   double GetALeft( void )  const { return tot1_; }
   double GetADown( void )  const { return tot2_; }
   double GetARight( void ) const { return tot2_; }
   
   double MeanTime( void )  const { return 0.5*(t1_+t2_); }
   double CMeanTime( void ) const { return 0.5*(ct1_+ct2_); }
   double DeltaE( void ) const;
   
   int DetectorId( void ) const { return raw_->DetectorId(); }
   int PlaneId( void ) const { return raw_->LayerId(); }
   int SegmentId( void ) const { return raw_->SegmentId(); }
   
   virtual bool ReCalc( bool applyRecursively=false )
      { return calculate(); }
   
};


#endif
