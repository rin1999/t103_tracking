/*
  ConfMan.cc

  2012/5 K.Shirotori
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "TrGeomMan.hh"
#include "FieldMan.hh"

const std::string defDCGeomFile="DCgeom.param";
// const std::string defFieldMapFile="FieldMap.dat";
const std::string defFieldMapFile="";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),
    anaMode_(false), 
    TrGeomFileName_(defDCGeomFile),
    TrGeomManager_(0),
    FieldMapFileName_(defFieldMapFile),
    SpecFieldNMR_(1.0), SpecFieldCalc_(1.0),
    BeamResol_(0.0),
    bSSDResol_(0.0), sSSD1Resol_(0.0), sSSD2Resol_(0.0), 
    BFTResol_(0.0), SFTResol_(0.0), AFTResol_(0.0),
    IT1Resol_(0.0), IT2Resol_(0.0), ST1Resol_(0.0), ST2Resol_(0.0), 
    T0Resol_(0.0), TofResol_(0.0), 
    ITofResol1_(0.0), ITofResol2_(0.0), PADResol_(0.0)
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
  if(TrGeomManager_) delete TrGeomManager_;
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
      //Analysis mode
      if( sscanf(buf,"ANAMODE: %d", &id )==1 ){
	if(id) anaMode_=true;
	else   anaMode_=false;
      }
      // Geometry
      else if( sscanf(buf,"TRGEO: %s",buf1)==1 ){
	TrGeomFileName_=buf1;
      }
      //Field Map
      else if( sscanf(buf,"FLDMAP: %s",buf1)==1 ){
	FieldMapFileName_=buf1;
      }
      else if( sscanf(buf,"FLDNMR: %lf", &val1)==1 )
	SpecFieldNMR_=val1;
      else if( sscanf(buf,"FLDCALC: %lf", &val1)==1 )
	SpecFieldCalc_=val1;
      //Resolution
      else if( sscanf(buf,"BEAMRESOL: %lf", &val1)==1 )
	BeamResol_=val1;
      else if( sscanf(buf,"bSSDRESOL: %lf", &val1)==1 )
	bSSDResol_=val1;
      else if( sscanf(buf,"sSSD1RESOL: %lf", &val1)==1 )
	sSSD1Resol_=val1;
      else if( sscanf(buf,"sSSD2RESOL: %lf", &val1)==1 )
	sSSD2Resol_=val1;
      else if( sscanf(buf,"BFTRESOL: %lf", &val1)==1 )
	BFTResol_=val1;
      else if( sscanf(buf,"SFTRESOL: %lf", &val1)==1 )
	SFTResol_=val1;
      else if( sscanf(buf,"AFTRESOL: %lf", &val1)==1 )
	AFTResol_=val1;
      else if( sscanf(buf,"IT1RESOL: %lf", &val1)==1 )
	IT1Resol_=val1;
      else if( sscanf(buf,"IT2RESOL: %lf", &val1)==1 )
	IT2Resol_=val1;
      else if( sscanf(buf,"ST1RESOL: %lf", &val1)==1 )
	ST1Resol_=val1;
      else if( sscanf(buf,"ST2RESOL: %lf", &val1)==1 )
	ST2Resol_=val1;
      else if( sscanf(buf,"T0RESOL: %lf", &val1)==1 )
	T0Resol_=val1;
      else if( sscanf(buf,"TOFRESOL: %lf", &val1)==1 )
	TofResol_=val1;
      else if( sscanf(buf,"ITOFRESOL1: %lf", &val1)==1 )
	ITofResol1_=val1;
      else if( sscanf(buf,"ITOFRESOL2: %lf", &val1)==1 )
	ITofResol2_=val1;
      else if( sscanf(buf,"PADRESOL: %lf", &val1)==1 )
	PADResol_=val1;

      else {
	std::cout << funcname << ": un-recognized record\n"
		  << buf << std::endl;
      }
    }
  }

  fclose(fp);
  InitializeParameterFiles();
  InitializeHistograms();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  TrGeomManager_ = & TrGeomMan::GetInstance();
  TrGeomManager_->Initialize(TrGeomFileName_);

  if( FieldMapFileName_!="" )
    FieldMan::GetInstance().Initialize( FieldMapFileName_ );

  return true;
}
