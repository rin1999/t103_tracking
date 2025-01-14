/*
  RawData.hh

  2018/10  K.Shirotori
*/

#ifndef RawData_h
#define RawData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DetectorID.hh"
#include <vector>

class Decoder;

class TrRawHit;
class HodoRawHit;
class HRTdcRawHit;

typedef std::vector<TrRawHit*>   TrRHitContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<HRTdcRawHit*> HRTdcRHitContainer;

typedef std::vector <double> DoubleVec;
typedef std::vector <unsigned int> uIntVec;
typedef std::vector <int> IntVec;

class RawData
{

private:  
  TrRHitContainer SFTRHC[PlMaxSFT+1];
  HodoRHitContainer T0RHC;
  HRTdcRHitContainer HRTRHC;

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeRawHits( Decoder& gDec );

  unsigned int getBigEndian32(const char* b);
  unsigned int Decode32bitWord(unsigned int word32bit);
  bool isAdcHg(unsigned int data);
  bool isAdcLg(unsigned int data);
  bool isTdcLeading(unsigned int data);
  bool isTdcTrailing(unsigned int data);
  bool isScaler(unsigned int data);

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddTrRHit( TrRHitContainer &cont, 
		  int Layer, int Fiber, int UorD, 
		  int AdcH, int AdcL, 
		  IntVec lTdc, IntVec tTdc );

  bool AddHodoRHit( HodoRHitContainer& cont,
		    int DetId, int Layer, int Seg, int UorD,
		    int EventNum,
		    DoubleVec Waveform,
		    DoubleVec Adc, DoubleVec Amp, 
		    DoubleVec Bl, DoubleVec PeakX,
		    uIntVec Tdc, uIntVec Dt, 
		    uIntVec Tdc_2nd, uIntVec Dt_2nd,
		    uIntVec Width, 
		    unsigned int L1_Tdc, unsigned int L1_Tdc1 
		    );

  bool AddHRTdcRHit( HRTdcRHitContainer& cont,
		     int DetId, int Layer, int Seg, int UorD,
		     int EventNum,
		     uIntVec lTdc, uIntVec tTdc
		     );

public: 
  const TrRHitContainer & GetSFTRHC( int layer ) const;
  const HodoRHitContainer& GetT0RHC() const;
  const HRTdcRHitContainer& GetHRTRHC() const;
};

#endif

