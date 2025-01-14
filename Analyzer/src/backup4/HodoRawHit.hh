/*
  HodoRawHit.hh
  
  2018/10  K.Shirotori
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <unsigned int> uIntVec;

class HodoRawHit
{
private:
  int DetId_, LayerId_, SegId_, UorD_;
  int EvNum_;
  DoubleVec Waveform_;
  DoubleVec Adc_, Amp_, Bl_, PeakX_;
  uIntVec Tdc_, Dt_, Tdc_2nd_, Dt_2nd_, Width_;
  unsigned int L1_Tdc_, L1_Tdc1_;
  
public:
  HodoRawHit( int detid, int layerid, int segid )
    : DetId_(detid), LayerId_(layerid), SegId_(segid),
      EvNum_(),
      UorD_(), Waveform_(),
      Adc_(), Amp_(), Bl_(), PeakX_(),
      Tdc_(), Dt_(), Tdc_2nd_(), Dt_2nd_(),
      Width_(),
      L1_Tdc_(), L1_Tdc1_()
  {};
  ~HodoRawHit() {};

public:
  int DetectorId( void ) const { return DetId_; };
  int LayerId( void ) const { return LayerId_; };
  int SegmentId( void ) const { return SegId_; };

  void SetEvNum( int evnum ) { EvNum_=evnum; }
  int GetEvNum( void ) const { return EvNum_; }

  void SetUorD( int uord ) { UorD_=uord; }
  int GetUorD( void ) const { return UorD_; }

  void SetWaveform( DoubleVec waveform ) { Waveform_=waveform; }

  void SetAdc( DoubleVec adc ) { Adc_=adc; }
  void SetAmp( DoubleVec amp ) { Amp_=amp; }
  void SetBl( DoubleVec bl ) { Bl_=bl; }
  void SetPeakX( DoubleVec peakx ) { PeakX_=peakx; }

  void SetTdc( uIntVec tdc ) { Tdc_=tdc; }
  void SetDt( uIntVec dt ) { Dt_=dt; }
  void SetTdc_2nd( uIntVec tdc2 ) { Tdc_2nd_=tdc2; }
  void SetDt_2nd( uIntVec dt2 ) { Dt_2nd_=dt2; }
  void SetWidth( uIntVec width ) { Width_=width; }

  void SetL1_Tdc( unsigned int l1_tdc ) { L1_Tdc_ =l1_tdc; }
  void SetL1_Tdc1( unsigned int l1_tdc1 ) { L1_Tdc1_=l1_tdc1; }

  double GetWaveform( int nh ) const { return Waveform_[nh]; };
  int  GetSize_Wf( void ) const { return Waveform_.size(); };

  double GetAdc( int nh ) const { return Adc_[nh]; };
  double GetAmp( int nh ) const { return Amp_[nh]; };
  double GetBl( int nh ) const { return Bl_[nh]; };
  double GetPeakX( int nh ) const { return PeakX_[nh]; };
  int  GetSize_Amp( void ) const { return Amp_.size(); };

  double GetTdc( int nh ) const { return Tdc_[nh]; };
  double GetDt( int nh ) const { return Dt_[nh]; };
  double GetTdc_2nd( int nh ) const { return Tdc_2nd_[nh]; };
  double GetDt_2nd( int nh ) const { return Dt_2nd_[nh]; };
  double GetWidth( int nh ) const { return Width_[nh]; };

  int  GetSize_Tdc( void ) const { return Tdc_.size(); };
  int  GetSize_Dt( void ) const { return Dt_.size(); };
  int  GetSize_Tdc_2nd( void ) const { return Tdc_2nd_.size(); };
  int  GetSize_Dt_2nd( void ) const { return Dt_2nd_.size(); };
  int  GetSize_Width( void ) const { return Width_.size(); };

  unsigned int GetL1_Tdc( void ) const { return L1_Tdc_; };
  unsigned int GetL1_Tdc1( void ) const { return L1_Tdc1_; };
};
#endif
