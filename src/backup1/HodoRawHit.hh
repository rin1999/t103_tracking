/*
  HodoRawHit.hh
  
  2016/2  K.Shirotori
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>
#include <vector>
#include <cfloat>

typedef std::vector <double> DoubleVec;
typedef std::vector <int> IntVec;

class HodoRawHit
{
private:
  int    DetId_, SegId_;
  DoubleVec TDC0_t_, TDC0_tot_;
  DoubleVec TDC1_t_, TDC1_tot_;
  DoubleVec ADC0_t_, ADC0_hgt_;
  DoubleVec ADC1_t_, ADC1_hgt_;

public:
  HodoRawHit( int detid, int segid )
    : DetId_(detid), SegId_(segid),
      TDC0_t_(0.0), TDC0_tot_(0.0),
      TDC1_t_(0.0), TDC1_tot_(0.0),
      ADC0_t_(0.0), ADC0_hgt_(0.0),
      ADC1_t_(0.0), ADC1_hgt_(0.0)
  {};
  ~HodoRawHit() {};

public:
  void SetTdc0Time( double t ) { TDC0_t_.push_back(t); }
  void SetTdc0Tot( double tot ) { TDC0_tot_.push_back(tot); }
  void SetTdc1Time( double t ) { TDC1_t_.push_back(t); }
  void SetTdc1Tot( double tot ) { TDC1_tot_.push_back(tot); }
  void SetAdc0Time( double t ) { ADC0_t_.push_back(t); }
  void SetAdc0Hgt( double hgt ) { ADC0_hgt_.push_back(hgt); }
  void SetAdc1Time( double t ) { ADC1_t_.push_back(t); }
  void SetAdc1Hgt( double hgt ) { ADC1_hgt_.push_back(hgt); }

  int DetectorId( void ) const { return DetId_; };
  int SegmentId( void ) const { return SegId_; };

  int    GetSize( void ) const { return TDC0_t_.size(); };
  double GetTDC0Time( int nh ) const { return TDC0_t_[nh]; };
  double GetTDC0Tot( int nh ) const { return TDC0_tot_[nh]; };
  double GetTDC1Time( int nh ) const { return TDC1_t_[nh]; };
  double GetTDC1Tot( int nh ) const { return TDC1_tot_[nh]; };
  double GetADC0Time( int nh ) const { return ADC0_t_[nh]; };
  double GetADC0Hgt( int nh ) const { return ADC0_hgt_[nh]; };
  double GetADC1Time( int nh ) const { return ADC1_t_[nh]; };
  double GetADC1Hgt( int nh ) const { return ADC1_hgt_[nh]; };
};
#endif
