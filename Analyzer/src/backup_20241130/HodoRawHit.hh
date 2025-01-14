/*
  HodoRawHit.hh
  
  2024/04  K.Shirotori
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class HodoRawHit
{

private:
   int DetId_, LayId_, SegId_;
   DoubleVec lTdc1_, Tot1_;
   DoubleVec lTdc2_, Tot2_;
   
public:
   HodoRawHit( int detid, int layid, int segid )
      : DetId_(detid),LayId_(layid), SegId_(segid),
        lTdc1_(0.0), Tot1_(0.0), lTdc2_(0.0), Tot2_(0.0)
      {};
   ~HodoRawHit() {};

public:
   int DetectorId( void ) const { return DetId_; };
   int LayerId( void ) const { return LayId_; };
   int SegmentId( void ) const { return SegId_; };
   
   void SetlTdc1( DoubleVec ltdc1 ) { lTdc1_=ltdc1; }
   void SetTot1( DoubleVec tot1 ) { Tot1_=tot1; }
   void SetlTdc2( DoubleVec ltdc2 ) { lTdc2_=ltdc2; }
   void SetTot2( DoubleVec tot2 ) { Tot2_=tot2; }

   double GetlTdc1( int nh ) const { return lTdc1_[nh]; }
   double GetTot1( int nh ) const { return Tot1_[nh]; }
   double GetlTdc2( int nh ) const { return lTdc2_[nh]; }
   double GetTot2( int nh ) const { return Tot2_[nh]; }
   
   int  GetSize_lTdc1( void ) const { return lTdc1_.size(); };
   int  GetSize_Tot1( void ) const { return Tot1_.size(); };
   int  GetSize_lTdc2( void ) const { return lTdc2_.size(); };
   int  GetSize_Tot2( void ) const { return Tot2_.size(); };

};
#endif
