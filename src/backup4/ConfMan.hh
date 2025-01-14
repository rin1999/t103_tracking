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
class TrTdcCalibMan;
class TrPHCMan;

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
  //TDC calibration                                                             
  std::string TrTdcCalibFileName_;
  TrTdcCalibMan *TrTdcCalibManager_;
  //Sleiwing correction (TOT)                                                   
  std::string TrPHCFileName_;
  TrPHCMan *TrPHCManager_;

  //TDC data
  bool HasTDC_;

  //Fiber analysis settings
  bool FiberTree_;
  int TRangeLow_, TRangeHigh_;
  int MinLayer_;

  //T0 analysis settings
  bool T0SimpleAna_;
  int T0RangeLow_, T0RangeHigh_;

  //DRS4 analysis settings
  int PeakStart_, PeakEnd_, BaseParam_;
  int RangeLow_, RangeHigh_, PileLow_, PileHigh_;

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
  //Tdc calib files                                                             
  TrTdcCalibMan *GetTrTdcCalibManager( void ) { return TrTdcCalibManager_; }
  //Slewing correction files                                                    
  TrPHCMan *GetTrPHCManager( void ) { return TrPHCManager_; }

  //TDC data
  bool HasTDC( void ) const { return HasTDC_; }
  //Fiber analysis settings
  bool FiberTree( void ) const { return FiberTree_; }
  int TRangeLow( void ) const { return TRangeLow_; }
  int TRangeHigh( void ) const { return TRangeHigh_; }
  int MinLayer( void ) const { return MinLayer_; }
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

private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

public:
  bool InitializeParameterFiles();
  bool InitializeHistograms();
};

#endif
