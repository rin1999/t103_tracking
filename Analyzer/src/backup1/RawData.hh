/*
  RawData.hh

  2016/2  K.Shirotori
*/

#ifndef RawData_h
#define RawData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DetectorID.hh"
#include <vector>

class PrimInfo;
class HodoRawHit;
class TrRawHit;

class s_BeamRawHit;
class s_ScatRawHit;

class TFile;

typedef std::vector<PrimInfo*>   PrimInfoContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<TrRawHit*>   TrRHitContainer;

typedef std::vector<s_BeamRawHit*> s_BeamRHitContainer;
typedef std::vector<s_ScatRawHit*> s_ScatRHitContainer;

class RawData
{

private:
  TrRHitContainer bSSDRHC[NumOfLayersbSSD+1];
  TrRHitContainer sSSD1RHC[NumOfLayerssSSD1+1];
  TrRHitContainer sSSD2RHC[NumOfLayerssSSD2+1];

  HodoRHitContainer T0RHC;
  HodoRHitContainer RPCRHC;

  s_BeamRHitContainer s_BeamRHC;
  s_ScatRHitContainer s_ScatRHC;

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeRawHits( TFile*, int );

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  //Full tracking
  bool AddTrRHit( TrRHitContainer &cont, 
		  int Layer, int Wire, 
		  double PosX, double PosY, double DL );
  
  bool AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int Seg,
			   double TDC0_t, double TDC0_tot,
			   double TDC1_t, double TDC1_tot,
			   double ADC0_t, double ADC0_hgt,
			   double ADC1_t, double ADC1_hgt);

  //Simple tracking
  bool AddsBTrRHit( s_BeamRHitContainer& cont, 
		    int TrackID, int Type,
		    int Layer, int Wire, 
		    double PosX, double PosY, double DL );
  
  bool AddsBHodoRHit( s_BeamRHitContainer& cont,
		      int TrackID, int Type,
		      int DetId, int Seg,
		      double TDC0_t, double TDC0_tot,
		      double TDC1_t, double TDC1_tot,
		      double ADC0_t, double ADC0_hgt,
		      double ADC1_t, double ADC1_hgt );

  bool AddsSTrRHit( s_ScatRHitContainer& cont, 
		    int TrackID, int Type, 
		    int Layer, int Wire, 
		    double PosX, double PosY, double DL );

  bool AddsSHodoRHit( s_ScatRHitContainer& cont,
		      int TrackID, int Type,
		      int DetId, int Seg,
		      double TDC0_t, double TDC0_tot,
		      double TDC1_t, double TDC1_tot,
		      double ADC0_t, double ADC0_hgt,
		      double ADC1_t, double ADC1_hgt );

public:
  const HodoRHitContainer& GetT0RHC() const;
  const HodoRHitContainer& GetRPCRHC() const;

  const TrRHitContainer & GetbSSDRHC( int layer ) const;
  const TrRHitContainer & GetsSSD1RHC( int layer ) const;
  const TrRHitContainer & GetsSSD2RHC( int layer ) const;

  const s_BeamRHitContainer & GetsBeamRHC() const;
  const s_ScatRHitContainer & GetsScatRHC() const;

};

#endif

