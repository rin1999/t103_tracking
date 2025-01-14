/*
  RawData.hh

  2024/04  K.Shirotori
*/

#ifndef RawData_h
#define RawData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DetectorInfo.hh"
#include <vector>

class HodoRawHit;
class DCRawHit;
class TrRawHit;

class TFile;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;
typedef std::vector<TrRawHit*>   TrRHitContainer;

typedef std::vector <double> DoubleVec;
typedef std::vector <int> IntVec;

class RawData
{

private:
   DCRHitContainer BDCRHC[NumOfLayersBDC+1];
   DCRHitContainer KLDCRHC[NumOfLayersKLDC+1];

   TrRHitContainer BFTRHC[NumOfLayersBFT+1];
   TrRHitContainer SFTRHC[NumOfLayersSFT+1];
   
   HodoRHitContainer UTOFRHC;
   HodoRHitContainer DTOFRHC;
   HodoRHitContainer LTOFRHC;
   HodoRHitContainer T0RHC;
   HodoRHitContainer T0rRHC;
   HodoRHitContainer BrefRHC;
   HodoRHitContainer T1RHC;
   HodoRHitContainer BHTRHC;
   
public:
   RawData();
   ~RawData();
   
   void clearAll();
   bool DecodeRawHits( TFile*, int );
   
private:
   RawData(const RawData&);
   RawData& operator=(const RawData&);
   
   bool AddDCRHit( DCRHitContainer &cont, 
                   int Layer, int Wire, DoubleVec lTdc, DoubleVec Tot );
   
   bool AddTrRHit( TrRHitContainer &cont, 
                   int Layer, int Fiber, DoubleVec lTdc, DoubleVec Tot );
   
   bool AddHodoRHit( HodoRHitContainer& cont,
                     int DetId, int Layer, int Seg,
                     DoubleVec lTdc1_, DoubleVec Tot1_,
                     DoubleVec lTdc2_, DoubleVec Tot2_ );
      
public:
   const DCRHitContainer & GetBDCRHC( int layer ) const;
   const DCRHitContainer & GetKLDCRHC( int layer ) const;

   const TrRHitContainer & GetBFTRHC( int layer ) const;
   const TrRHitContainer & GetSFTRHC( int layer ) const;
      
   const HodoRHitContainer& GetUTOFRHC() const;
   const HodoRHitContainer& GetDTOFRHC() const;
   const HodoRHitContainer& GetLTOFRHC() const;
   const HodoRHitContainer& GetT0RHC() const;
   const HodoRHitContainer& GetT0rRHC() const;
   const HodoRHitContainer& GetBrefRHC() const;
   const HodoRHitContainer& GetT1RHC() const;
   const HodoRHitContainer& GetBHTRHC() const;
};

#endif

