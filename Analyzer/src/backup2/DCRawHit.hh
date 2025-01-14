/*
  DCRawHit.hh
  
  2018/12  K.Shirotori
*/

#ifndef DCRawHit_h 
#define DCRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class TrLTrackHit;

class DCRawHit
{
private:
  int LayerId_, WireId_;
  IntVec lTdc_, tTdc_;
  
public:
  DCRawHit( int layer, int wire )
    : LayerId_(layer), WireId_(wire),
      lTdc_(0.0), tTdc_(0.0)
  {};
  ~DCRawHit() {};

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetWire( int wire ) { WireId_=wire; }
  void SetlTdc( IntVec rtdc ) { lTdc_=rtdc; }
  void SettTdc( IntVec ttdc ) { tTdc_=ttdc; }

  int LayerId( void ) const { return LayerId_; }
  int WireId( void ) const { return WireId_; }
  double GetlTdc( int nh ) const { return lTdc_[nh]; }
  double GettTdc( int nh ) const { return tTdc_[nh]; }

  int GetSize_lTdc( void ) const { return lTdc_.size(); }
  int GetSize_tTdc( void ) const { return tTdc_.size(); }
};
#endif
