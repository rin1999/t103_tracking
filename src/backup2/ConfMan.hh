/*
  ConfMan.hh

  2018/10 K.Shirotori
*/

#ifndef Confman_h
#define Confman_h 1

#include <string>

class VEvent;
class CMapMan;
class TrGeomMan;
class DCGeomMan;
class DCTdcCalibMan;
class DCDriftParamMan;

class ConfMan
{

private:
  //Conf
  std::string ConfFileName_;
  static ConfMan *confManager_;

  //Counter map
  std::string CMapFileName_;
  CMapMan *CMapManager_;

  //Geometry files  
  std::string TrGeomFileName_;
  TrGeomMan *TrGeomManager_;
  DCGeomMan *DCGeomManager_;
  std::string      DCTdcCalibFileName_;
  DCTdcCalibMan   *DCTdcCalibManager_;
  std::string      DCDriftParamFileName_;
  DCDriftParamMan *DCDriftParamManager_;

  //TDC data
  bool HasTDC_;

  //T0 analysis settings
  bool T0SimpleAna_;
  int T0RangeLow_, T0RangeHigh_;

  //DRS4 analysis settings
  int PeakStart_, PeakEnd_, BaseParam_;
  int RangeLow_, RangeHigh_, PileLow_, PileHigh_;

  //DC analysis settings
  bool DCTree_;
  int DCTRangeLow_, DCTRangeHigh_;
  int DCWidthCutL_, DCWidthCutH_;
  double DTCut1_, DTCut2_, DTCut3_, DTCut4_;
  double DLCut1_, DLCut2_, DLCut3_, DLCut4_;

public:
  ConfMan( const std::string & filename );
  ~ConfMan();

  static ConfMan *GetConfManager( void ) { return confManager_; }
  bool Initialize( void );
  void SetFileName( const std::string & filename ) { ConfFileName_=filename; }
  VEvent* EventAllocator();

  //Counter map
  CMapMan *GetCMapManager( void )   { return CMapManager_; }
 //Geometry files  
  TrGeomMan *GetTrGeomManager( void ) { return TrGeomManager_; }
  DCGeomMan *GetDCGeomManager( void ) { return DCGeomManager_; }
  DCTdcCalibMan   *GetDCTdcCalibManager( void ) { return DCTdcCalibManager_; }
  DCDriftParamMan *GetDCDriftParamManager( void ) { return DCDriftParamManager_; } 

  //TDC data
  bool HasTDC( void ) const { return HasTDC_; }

  //T0 analysis settings
  bool T0SimpleAna( void ) const { return T0SimpleAna_; }
  int T0RangeLow( void ) const { return T0RangeLow_; }
  int T0RangeHigh( void ) const { return T0RangeHigh_; }
  //DRS4 analysis settings
  int PeakStart( void ) const { return PeakStart_; }
  int PeakEnd( void ) const { return PeakEnd_; }
  int BaseParam( void ) const { return BaseParam_; }
  int RangeLow( void ) const { return RangeLow_; }
  int RangeHigh( void ) const { return RangeHigh_; }
  int PileLow( void ) const { return PileLow_; }
  int PileHigh( void ) const { return PileHigh_; }
  //DC analysis settings
  bool DCTree( void ) const { return DCTree_; }
  int DCTRangeLow( void ) const { return DCTRangeLow_; }
  int DCTRangeHigh( void ) const { return DCTRangeHigh_; }
  int DCWidthCutLow( void ) const { return DCWidthCutL_; }
  int DCWidthCutHigh( void ) const { return DCWidthCutH_; }
  double DTCut1( void ) const { return DTCut1_; }
  double DTCut2( void ) const { return DTCut2_; }
  double DTCut3( void ) const { return DTCut3_; }
  double DTCut4( void ) const { return DTCut4_; }
  double DLCut1( void ) const { return DLCut1_; }
  double DLCut2( void ) const { return DLCut2_; }
  double DLCut3( void ) const { return DLCut3_; }
  double DLCut4( void ) const { return DLCut4_; }

private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

public:
  bool InitializeParameterFiles();
  bool InitializeHistograms();
};

#endif
