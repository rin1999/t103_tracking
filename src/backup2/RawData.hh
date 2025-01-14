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

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;

typedef std::vector <double> DoubleVec;
typedef std::vector <unsigned int> uIntVec;
typedef std::vector <int> IntVec;

class RawData
{

private:
  HodoRHitContainer T0RHC;
  DCRHitContainer DCRHC[PlMaxDC+1];

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeRawHits( Decoder& gDec );

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddHodoRHit( HodoRHitContainer& cont,
		    int DetId, int Layer, int Seg, int UorD,
		    uIntVec lTdc, uIntVec tTdc, uIntVec Width );

  bool AddDCRHit( DCRHitContainer &cont, 
		  int Layer, int Wire, IntVec lTdc, IntVec tTdc );

public:
  const HodoRHitContainer& GetT0RHC() const;
  const DCRHitContainer & GetDCRHC( int layer ) const;
};

#endif

