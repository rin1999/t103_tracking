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
  int  Adc_;
  IntVec lTdc_, tTdc_;
  
public:
  TrRawHit( int layer, int fiber, int uord )
    : LayerId_(layer), FiberId_(fiber), UorD_(uord),
      Adc_(0), lTdc_(0.0), tTdc_(0.0)
  {};
  ~TrRawHit() {};

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetFiber( int fiber ) { FiberId_=fiber; }
  void SetUorD( int ud ) { UorD_=ud; }
  void SetAdc( int adc ) { Adc_=adc; }
  void SetlTdc( IntVec rtdc ) { lTdc_=rtdc; }
  void SettTdc( IntVec ttdc ) { tTdc_=ttdc; }

  int LayerId( void ) const { return LayerId_; }
  int FiberId( void ) const { return FiberId_; }
  int GetUorD( void ) const { return UorD_; }
  int GetAdc( void ) const { return Adc_; }
  double GetlTdc( int nh ) const { return lTdc_[nh]; }
  double GettTdc( int nh ) const { return tTdc_[nh]; }

  int GetSize_lTdc( void ) const { return lTdc_.size(); }
  int GetSize_tTdc( void ) const { return tTdc_.size(); }
};
#endif
