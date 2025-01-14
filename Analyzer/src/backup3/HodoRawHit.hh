/*
 HodoRawHit.hh
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>

class HodoRawHit
{

private:
  int DetId_, PlId_, SegId_;
  int Adc1_, Adc2_;
  int Tdc1_, Tdc2_;
  int NhitsTdc_;

public:
  HodoRawHit( int detid, int plid, int segid )
    : DetId_(detid), PlId_(plid), SegId_(segid),
      Adc1_(-1), Adc2_(-1), Tdc1_(-1), Tdc2_(-1), NhitsTdc_(0)
  {};

  ~HodoRawHit() {};

public:
  void SetAdc1(int adc) { Adc1_=adc; };
  void SetAdc2(int adc) { Adc2_=adc; };
  void SetTdc1(int tdc) {
    if( Tdc1_==-1 ) { Tdc1_=tdc; ++NhitsTdc_; }
    else    Tdc1_=tdc;
  };
  void SetTdc2(int tdc) {
    if( Tdc2_==-1 ) { Tdc2_=tdc; ++NhitsTdc_; }
    else    Tdc2_=tdc;
  };

  void SetAdcUp( int adc) { SetAdc1(adc); };
  void SetAdcLeft( int adc) { SetAdc1(adc); };
  void SetAdcDown( int adc) { SetAdc2(adc); };
  void SetAdcRight( int adc) { SetAdc2(adc); };
  void SetTdcUp( int tdc) { SetTdc1(tdc); };
  void SetTdcLeft( int tdc) { SetTdc1(tdc); };
  void SetTdcDown( int tdc) { SetTdc2(tdc); };
  void SetTdcRight( int tdc) { SetTdc2(tdc); };

  int DetectorId( void ) const { return DetId_; };
  int PlaneId( void ) const { return PlId_; };
  int SegmentId( void ) const { return SegId_; };

  int GetAdc1( void ) const { return Adc1_; };
  int GetAdc2( void ) const { return Adc2_; };
  int GetTdc1( void ) const { return Tdc1_; };
  int GetTdc2( void ) const { return Tdc2_; };

  int GetNumOfTdcHits( void ) const { return NhitsTdc_; };

  int GetAdcUp( void ) const { return GetAdc1(); };
  int GetAdcLeft( void ) const { return GetAdc1(); };
  int GetAdcDown( void ) const { return GetAdc2(); };
  int GetAdcRight( void ) const { return GetAdc2(); };
  int GetTdcUp( void ) const { return GetTdc1(); };
  int GetTdcLeft( void ) const { return GetTdc1(); };
  int GetTdcDown( void ) const { return GetTdc2(); };
  int GetTdcRight( void ) const { return GetTdc2(); };
};
#endif
