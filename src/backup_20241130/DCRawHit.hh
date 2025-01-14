/*
  DCRawHit.hh
  
  2024/04  K.Shirotori
*/

#ifndef DCRawHit_h 
#define DCRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class DCLTrackHit;

class DCRawHit
{
private:
  int LayerId_, WireId_;
  DoubleVec lTdc_, Tot_;
  
public:
  DCRawHit( int layerid, int wireid )
    : LayerId_(layerid), WireId_(wireid),
      lTdc_(0.0), Tot_(0.0)
  {};
  ~DCRawHit() {};

public:   
  void SetLayer( int layerid ) { LayerId_=layerid; }
  void SetWire( int wireid ) { WireId_=wireid; }
  void SetlTdc( DoubleVec ltdc ) { lTdc_=ltdc; }
  void SetTot( DoubleVec tot ) { Tot_=tot; }

  int LayerId( void ) const { return LayerId_; }
  int WireId( void ) const { return WireId_; }
  double GetlTdc( int nh ) const { return lTdc_[nh]; }
  double GetTot( int nh ) const { return Tot_[nh]; }

  int GetSize_lTdc( void ) const { return lTdc_.size(); }
  int GetSize_Tot( void ) const { return Tot_.size(); }
};
#endif
