/*
  ConfMan.cc

  2018/10 K.Shirotori
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "CMapMan.hh"
#include "TrGeomMan.hh"
#include "TrTdcCalibMan.hh"
#include "TrPHCMan.hh"

const std::string defCMapFile="CMap.param";
const std::string defTrGeomFile="TrGeom.param";
const std::string defTrTdcCalibFile="TrTdcCalib.param";
const std::string defTrPHCFile="TrPHC.param";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),
    CMapManager_(0), 
    TrGeomFileName_(defTrGeomFile), TrGeomManager_(0),
    TrTdcCalibFileName_(defTrTdcCalibFile), TrTdcCalibManager_(0),
    TrPHCFileName_(defTrPHCFile), TrPHCManager_(0),
    HasTDC_(0),
    FiberTree_(0), TRangeLow_(0), TRangeHigh_(0), MinLayer_(12),
    T0SimpleAna_(0), T0RangeLow_(0), T0RangeHigh_(0),
    PeakStart_(0), PeakEnd_(0), BaseParam_(0), 
    RangeLow_(0), RangeHigh_(0), PileLow_(0), PileHigh_(0)
{
  static const std::string funcname="[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constring twice" << std::endl;
    std::exit(-1);
  }
  confManager_=this;
}

ConfMan::~ConfMan( )
{
  confManager_=0;
}

const int BufSize=300;

bool ConfMan::Initialize()
{
  static const std::string funcname="[ConfMan::Initialize]";

  char buf[BufSize], buf1[BufSize];
  double val1, val2, val3;
  int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){ 
      //Counter map
      if( sscanf(buf,"CMAP: %s",buf1)==1 ){
	CMapFileName_=buf1;
      }
      // Geometry                                                               
      else if( sscanf(buf,"TRGEO: %s",buf1)==1 ){
        TrGeomFileName_=buf1;
      }
      // Tdc Calibration                                                        
      else if( sscanf(buf,"TRTDC: %s",buf1)==1 ){
        TrTdcCalibFileName_=buf1;
      }
      // Slewing correction                                                     
      else if( sscanf(buf,"TRPHC: %s",buf1)==1 ){
        TrPHCFileName_=buf1;
      }
      //DAQ conditions
      else if( sscanf(buf,"HASTDC: %d", &id )==1 ){
	HasTDC_=id;
      }
      //Fiber
      else if( sscanf(buf,"FIBERTREE: %d", &id )==1 ){
        FiberTree_=id;
      }
      else if( sscanf(buf,"TRANGELOW: %d", &id )==1 ){
	TRangeLow_=id;
      }
      else if( sscanf(buf,"TRANGEHIGH: %d", &id )==1 ){
	TRangeHigh_=id;
      }
      else if( sscanf(buf,"MINLAYER: %d", &id )==1 ){
	MinLayer_=id;
      }
      //T0
      else if( sscanf(buf,"T0SIMPLEANA: %d", &id )==1 ){
	T0SimpleAna_=id;
      }
      else if( sscanf(buf,"T0RANGELOW: %d", &id )==1 ){
	T0RangeLow_=id;
      }
      else if( sscanf(buf,"T0RANGEHIGH: %d", &id )==1 ){
	T0RangeHigh_=id;
      }
      //DRS4
      else if( sscanf(buf,"PEAKSTART: %d", &id )==1 ){
	PeakStart_=id;
      }
      else if( sscanf(buf,"PEAKEND: %d", &id )==1 ){
	PeakEnd_=id;
      }
      else if( sscanf(buf,"BASEPARAM: %d", &id )==1 ){
	BaseParam_=id;
      }
      else if( sscanf(buf,"RANGELOW: %d", &id )==1 ){
	RangeLow_=id;
      }
      else if( sscanf(buf,"RANGEHIGH: %d", &id )==1 ){
	RangeHigh_=id;
      }
      else if( sscanf(buf,"PILELOW: %d", &id )==1 ){
	PileLow_=id;
      }
      else if( sscanf(buf,"PILEHIGH: %d", &id )==1 ){
	PileHigh_=id;
      }
      else {
	std::cout << funcname << ": un-recognized record\n"
		  << buf << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << "------------------" << ConfFileName_ << "------" << std::endl;
  std::cout << "Counter Map:      " << CMapFileName_ << std::endl;
  std::cout << "Tracker Geom. Param.:  " << TrGeomFileName_ << std::endl;
  std::cout << "Tracker TdcCalib Param.:  " << TrTdcCalibFileName_ << std::endl;
  std::cout << "Tracker PHC Param.:  " << TrPHCFileName_ << std::endl;
  std::cout << "DRS4 TDC decoding:  " << HasTDC_ << std::endl;
  std::cout << "***** Analysis condistions" << std::endl;
  std::cout << "Fiber tree construction:  " << FiberTree_ << std::endl;
  std::cout << "Fiber TDC range:  " << TRangeLow_ << " < TDC <" << TRangeHigh_ << std::endl;
  std::cout << "Fiber Minimum Layer hit:  " << MinLayer_ << std::endl;
  std::cout << "T0 simple analysis:  " << T0SimpleAna_ << std::endl;
  std::cout << "T0 TDC range:  " << T0RangeLow_ << " < TDC <" << T0RangeHigh_ << std::endl;
  std::cout << "***** DRS4 condistions" << std::endl;
  std::cout << "Peak Search:  " << PeakStart_ << " - " << PeakEnd_ << std::endl;
  std::cout << "Base line:  " << BaseParam_ << std::endl;
  std::cout << "Integral:  " << RangeLow_ << " - " << RangeHigh_ << std::endl;
  std::cout << "Pileup:  " << PileLow_ << " - " << PileHigh_ << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  InitializeParameterFiles();
  InitializeHistograms();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  if( CMapFileName_!="" )
    CMapManager_ = new CMapMan(CMapFileName_);
  if(CMapManager_) CMapManager_->Initialize();

  TrGeomManager_ = & TrGeomMan::GetInstance();
  if( TrGeomFileName_!="" )
    TrGeomManager_->Initialize(TrGeomFileName_);
  else
    TrGeomManager_->Initialize();

  if( TrTdcCalibFileName_!="" )
    TrTdcCalibManager_ = new TrTdcCalibMan(TrTdcCalibFileName_);
  if(TrTdcCalibManager_) TrTdcCalibManager_->Initialize();

  if( TrPHCFileName_!="" )
    TrPHCManager_ = new TrPHCMan(TrPHCFileName_);
  if(TrPHCManager_) TrPHCManager_->Initialize();

  return true;
}
