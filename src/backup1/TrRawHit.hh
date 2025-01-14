/*
  TrRawHit.hh
  
  2016/2  K.Shirotori
*/

#ifndef TrRawHit_h 
#define TrRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <bool> BoolVec;

class TrLTrackHit;

class TrRawHit
{
private:
  int LayerId_, WireId_;
  DoubleVec  PosX_, PosY_, DL_;

public:
  TrRawHit( int layer, int wire )
    : LayerId_(layer), WireId_(wire), 
      PosX_(0.0), PosY_(0.0), DL_(0.0)
  {};
  ~TrRawHit() {};

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetWire( int wire ) { WireId_=wire; }
  void SetPosX( double posx ) { PosX_.push_back(posx); }
  void SetPosY( double posy ) { PosY_.push_back(posy); }
  void SetDL( double dl ) { DL_.push_back(dl); }

  int LayerId( void ) const { return LayerId_; }
  int WireId( void ) const { return WireId_; }
  double GetPosX( int nh ) const { return PosX_[nh]; }
  double GetPosY( int nh ) const { return PosY_[nh]; }
  double GetDL( int nh ) const { return DL_[nh]; }
  int GetSize( void ) const { return DL_.size(); }

};
#endif
