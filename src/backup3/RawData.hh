/*
  RawData.hh
*/

#ifndef RawData_h
#define RawData_h

#include "DetectorID.hh"
#include <vector>

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector <DCRawHit *> DCRHitContainer;

class RawData
{

private:
  HodoRHitContainer GCRawHC;
  HodoRHitContainer BH1RawHC;
  HodoRHitContainer BH2RawHC;
  HodoRHitContainer BACRawHC;
  HodoRHitContainer TGTRawHC;

  HodoRHitContainer TOFRawHC;
  HodoRHitContainer LCRawHC;
  HodoRHitContainer ACRawHC[NumOfLayersAc+1];
  HodoRHitContainer SP0RawHC[NumOfLayersSP0+1];

  DCRHitContainer BcInRawHC[NumOfLayersBcIn+1],   BcOutRawHC[NumOfLayersBcOut+1];
  DCRHitContainer SdcInRawHC[NumOfLayersSdcIn+1], SdcOutRawHC[NumOfLayersSdcOut+1];

  HodoRHitContainer MiscRawHC;
  HodoRHitContainer MatrixRawHC;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeHits();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddHodoRawHit(HodoRHitContainer& cont,
 		     int DetId,
 		     int Plane,
 		     int Seg,
 		     int AorT,
 		     int UorD,
 		     int Data);

  bool AddDCRawHit( DCRHitContainer &cont, 
		    int Plane, 
		    int Wire, 
		    int Tdc,
		    int type=0); 

public:
  const HodoRHitContainer& GetGCRawHC() const;
  const HodoRHitContainer& GetBH1RawHC() const;
  const HodoRHitContainer& GetBH2RawHC() const;
  const HodoRHitContainer& GetBACRawHC() const;
  const HodoRHitContainer& GetTGTRawHC() const;

  const HodoRHitContainer& GetTOFRawHC() const;
  const HodoRHitContainer& GetLCRawHC() const;
  const HodoRHitContainer& GetACRawHC( int layer ) const;
  const HodoRHitContainer& GetSP0RawHC( int layer ) const;

  const DCRHitContainer & GetBcInRawHC( int layer ) const;
  const DCRHitContainer & GetBcOutRawHC( int layer ) const;
  const DCRHitContainer & GetSdcInRawHC( int layer ) const;
  const DCRHitContainer & GetSdcOutRawHC( int layer ) const;

  const HodoRHitContainer& GetMiscRawHC() const;
  const HodoRHitContainer& GetMatrixRawHC() const;
};

#endif
