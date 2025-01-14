/*
  TrRawHit.hh
  
  2024/10  K.Shirotori
*/

#ifndef TrRawHit_h 
#define TrRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class TrLTrackHit;

class TrRawHit
{
private:
  int LayerId_, FiberId_;
  DoubleVec lTdc_, Tot_;
  
public:
  TrRawHit( int layerid, int fiberid )
    : LayerId_(layerid), FiberId_(fiberid),
      lTdc_(0.0), Tot_(0.0)
  {};
  ~TrRawHit() {};

public:   
  void SetLayer( int layerid ) { LayerId_=layerid; }
  void SetFiber( int fiberid ) { FiberId_=fiberid; }
  void SetlTdc( DoubleVec ltdc ) { lTdc_=ltdc; }
  void SetTot( DoubleVec tot ) { Tot_=tot; }

  int LayerId( void ) const { return LayerId_; }
  int FiberId( void ) const { return FiberId_; }
  double GetlTdc( int nh ) const { return lTdc_[nh]; }
  double GetTot( int nh ) const { return Tot_[nh]; }

  int GetSize_lTdc( void ) const { return lTdc_.size(); }
  int GetSize_Tot( void ) const { return Tot_.size(); }
};
#endif
