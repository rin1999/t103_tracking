/*
  TrRawHit.hh
  
  2018/10  K.Shirotori
*/

#ifndef TrRawHit_h 
#define TrRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;
typedef std::vector <bool> BoolVec;

class TrLTrackHit;

class TrRawHit
{
private:
  int LayerId_, FiberId_, UorD_;
  int  AdcH_, AdcL_;
  IntVec lTdc_, tTdc_;
  
public:
  TrRawHit( int layer, int fiber, int uord )
    : LayerId_(layer), FiberId_(fiber), UorD_(uord),
      AdcH_(0), AdcL_(0), lTdc_(0.0), tTdc_(0.0)
  {};
  ~TrRawHit() {};

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetFiber( int fiber ) { FiberId_=fiber; }
  void SetUorD( int ud ) { UorD_=ud; }
  void SetAdcH( int adch ) { AdcH_=adch; }
  void SetAdcL( int adcl ) { AdcL_=adcl; }
  void SetlTdc( IntVec rtdc ) { lTdc_=rtdc; }
  void SettTdc( IntVec ttdc ) { tTdc_=ttdc; }

  int LayerId( void ) const { return LayerId_; }
  int FiberId( void ) const { return FiberId_; }
  int GetUorD( void ) const { return UorD_; }
  int GetAdcH( void ) const { return AdcH_; }
  int GetAdcL( void ) const { return AdcL_; }
  double GetlTdc( int nh ) const { return lTdc_[nh]; }
  double GettTdc( int nh ) const { return tTdc_[nh]; }

  int GetSize_lTdc( void ) const { return lTdc_.size(); }
  int GetSize_tTdc( void ) const { return tTdc_.size(); }
};
#endif
