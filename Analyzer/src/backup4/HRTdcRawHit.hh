/*
  HRTdcRawHit.hh
  
  2019/12  K.Shirotori
*/

#ifndef HRTdcRawHit_h 
#define HRTdcRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <unsigned int> uIntVec;

class HRTdcRawHit
{
private:
  int DetId_, LayerId_, SegId_, UorD_;
  int EvNum_;
  uIntVec lTdc_, tTdc_;
  
public:
  HRTdcRawHit( int detid, int layerid, int segid )
    : DetId_(detid), LayerId_(layerid), SegId_(segid),
      EvNum_(),
      UorD_(), 
      lTdc_(), tTdc_()
  {};
  ~HRTdcRawHit() {};
  
public:
  int DetectorId( void ) const { return DetId_; };
  int LayerId( void ) const { return LayerId_; };
  int SegmentId( void ) const { return SegId_; };

  void SetEvNum( int evnum ) { EvNum_=evnum; }
  int GetEvNum( void ) const { return EvNum_; }

  void SetUorD( int uord ) { UorD_=uord; }
  int GetUorD( void ) const { return UorD_; }

  void SetlTdc( uIntVec ltdc ) { lTdc_=ltdc; }
  void SettTdc( uIntVec ttdc ) { tTdc_=ttdc; }

  double GetlTdc( int nh ) const { return lTdc_[nh]; };
  double GettTdc( int nh ) const { return tTdc_[nh]; };

  int  GetSize_lTdc( void ) const { return lTdc_.size(); };
  int  GetSize_tTdc( void ) const { return tTdc_.size(); };
};
#endif
