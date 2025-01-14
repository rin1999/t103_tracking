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
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"

const std::string defCMapFile="CMap.param";
const std::string defDCGeomFile="DCgeom.param";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),
    CMapManager_(0), 
    TrGeomFileName_(defDCGeomFile), TrGeomManager_(0),
    DCTdcCalibManager_(0),
    DCDriftParamManager_(0),
    HasTDC_(0),
    T0SimpleAna_(0), T0RangeLow_(0), T0RangeHigh_(0),
    PeakStart_(0), PeakEnd_(0), BaseParam_(0), 
    RangeLow_(0), RangeHigh_(0), PileLow_(0), PileHigh_(0),
    DCTree_(0), DCTRangeLow_(0), DCTRangeHigh_(0), 
    DCWidthCutL_(0), DCWidthCutH_(0),
    DTCut1_(0), DTCut2_(0), DTCut3_(0), DTCut4_(0), 
    DLCut1_(0), DLCut2_(0), DLCut3_(0), DLCut4_(0)
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
      else if( sscanf(buf,"DCTDC: %s",buf1)==1 ){
	DCTdcCalibFileName_=buf1;
      }
      else if( sscanf(buf,"DCDRFT: %s",buf1)==1 ){
	DCDriftParamFileName_=buf1;
      }
      //DAQ conditions
      else if( sscanf(buf,"HASTDC: %d", &id )==1 ){
	HasTDC_=id;
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
      //DC
      else if( sscanf(buf,"DCTREE: %d", &id )==1 ){
	DCTree_=id;
      }
      else if( sscanf(buf,"DCTRANGELOW: %d", &id )==1 ){
	DCTRangeLow_=id;
      }
      else if( sscanf(buf,"DCTRANGEHIGH: %d", &id )==1 ){
	DCTRangeHigh_=id;
      }
      else if( sscanf(buf,"DCWIDTHCUTL: %d", &id )==1 ){
	DCWidthCutL_=id;
      }
      else if( sscanf(buf,"DCWIDTHCUTH: %d", &id )==1 ){
	DCWidthCutH_=id;
      }
      else if( sscanf(buf,"DTCut1: %lf", &val1 )==1 ){
	DTCut1_=val1;
      }
      else if( sscanf(buf,"DTCut2: %lf", &val1 )==1 ){
	DTCut2_=val1;
      }
      else if( sscanf(buf,"DTCut3: %lf", &val1 )==1 ){
	DTCut3_=val1;
      }
      else if( sscanf(buf,"DTCut4: %lf", &val1 )==1 ){
	DTCut4_=val1;
      }
      else if( sscanf(buf,"DLCut1: %lf", &val1 )==1 ){
	DLCut1_=val1;
      }
      else if( sscanf(buf,"DLCut2: %lf", &val1 )==1 ){
	DLCut2_=val1;
      }
      else if( sscanf(buf,"DLCut3: %lf", &val1 )==1 ){
	DLCut3_=val1;
      }
      else if( sscanf(buf,"DLCut4: %lf", &val1 )==1 ){
	DLCut4_=val1;
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
  std::cout << "DC TDC Calib. Param.:  " << DCTdcCalibFileName_ << std::endl;
  std::cout << "DC Drift Time Param.:  " << DCDriftParamFileName_ << std::endl;
  std::cout << "DRS4 TDC decoding:  " << HasTDC_ << std::endl;
  std::cout << "***** Analysis condistions" << std::endl;
  std::cout << "T0 simple analysis:  " << T0SimpleAna_ << std::endl;
  std::cout << "T0 TDC range:  " << T0RangeLow_ << " < TDC <" << T0RangeHigh_ << std::endl;
  std::cout << "***** DRS4 condistions" << std::endl;
  std::cout << "Peak Search:  " << PeakStart_ << " - " << PeakEnd_ << std::endl;
  std::cout << "Base line:  " << BaseParam_ << std::endl;
  std::cout << "Integral:  " << RangeLow_ << " - " << RangeHigh_ << std::endl;
  std::cout << "Pileup:  " << PileLow_ << " - " << PileHigh_ << std::endl;
  std::cout << "***** DC condistions" << std::endl;
  std::cout << "DC tree construction:  " << DCTree_ << std::endl;
  std::cout << "DC TDC range:  " << DCTRangeLow_ << " < TDC <" << DCTRangeHigh_ << std::endl;
  std::cout << "DC Width Cut:  " << DCWidthCutL_ << " < Width <" << DCWidthCutH_ << std::endl;
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

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( TrGeomFileName_!="" )
    DCGeomManager_->Initialize(TrGeomFileName_);
  else
    DCGeomManager_->Initialize();

  if( DCTdcCalibFileName_!="" )
    DCTdcCalibManager_ = new DCTdcCalibMan(DCTdcCalibFileName_);
  if(DCTdcCalibManager_) DCTdcCalibManager_->Initialize();

  if( DCDriftParamFileName_!="" )
    DCDriftParamManager_ = new DCDriftParamMan(DCDriftParamFileName_);
  if(DCDriftParamManager_) DCDriftParamManager_->Initialize();

  return true;
}
