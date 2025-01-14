/*
  ConfMan.cc
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>

#include "UnpackerManager.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "K18TransMatrix.hh"

#include "EvDisp.hh"

ConfMan * ConfMan::confManager_=0;

ConfMan::ConfMan()
  : m_ConfFileName(),
    m_UnpackerFileName(),
    m_DigitFileName(),
    m_CmapFileName(),
    DCGeomManager_(0),
    DCTdcCalibManager_(0),
    DCDriftParamManager_(0),
    K18Matrix_(),
    K18Momentum_(1.05),
    SKSFieldNMR_(2.2), SKSFieldCalc_(2.2),
    evDisp_(0), FlagEvDisp_(false),
    FlagMHTDC_(false)
{
}

ConfMan::ConfMan(const std::string& confFile)
  :m_ConfFileName(confFile),
   m_UnpackerFileName()
{
  if (confManager_)
    std::exit(-1);
    
  confManager_ = this;
}

ConfMan::~ConfMan()
{
}

// ConfMan&
// ConfMan::GetInstance()
// {
//   static ConfMan gConfMan;
//   return gConfMan;
// }

bool
ConfMan::Initialize()
{
  double val;
  std::ifstream f(m_ConfFileName.c_str());
  //m_ConfFileName = confFile;
  if (f.fail()){
    std::cerr << "#E " 
	      << std::endl;
    return false;
  }
  
  while (f.good()){
    std::string line;
    std::getline(f, line);
    if (line.empty())
      continue;
    std::istringstream input_line(line);
    std::istream_iterator<std::string> line_begin(input_line);
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> param_list(line_begin, line_end);
    if (param_list.size()!=2)
      continue;
    std::cout << param_list[0] << " " << param_list[1] << std::endl;
    //Conf
    if (param_list[0] == "UNPACKER:")
      m_UnpackerFileName = param_list[1];
    if (param_list[0] == "DIGIT:")
      m_DigitFileName = param_list[1];
    if (param_list[0] == "CMAP:")
      m_CmapFileName = param_list[1];
    //Hodo
    if (param_list[0] == "HDPRM:")
      HodoParamFileName_ = param_list[1];	
    if (param_list[0] == "HDPHC:")
      HodoPHCFileName_ = param_list[1];	
    //DC
    if (param_list[0] == "DCGEO:")
      DCGeomFileName_ = param_list[1];
    if (param_list[0] == "DCTDC:")
      DCTdcCalibFileName_ = param_list[1];
    if (param_list[0] == "DCDRFT:")
      DCDriftParamFileName_ = param_list[1];
    if (param_list[0] == "K18TM:")
      K18MatrixFileName_ = param_list[1];
    if (param_list[0] == "PK18:")
      K18Momentum_ = boost::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDMAP:")
      FieldMapFileName_ = param_list[1];
    if (param_list[0] == "FLDNMR:")
      SKSFieldNMR_ = boost::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDCALC:")
      SKSFieldCalc_ = boost::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "BH1FLT:")
      bh1FilterFileName_ = param_list[1];
//     if (param_list[0] == "D4FUNC:")
//       D4FuncFileName_ = param_list[1];

    // Event display 
    if( param_list[0] == "EVDISP:" )
      if(param_list[1]=="1") FlagEvDisp_=true;
      else         FlagEvDisp_=false;
      
    //MHTDC event suppression flag
    if( param_list[0] == "MHTDC:" )
      if(param_list[1]=="1") FlagMHTDC_=true;
      else         FlagMHTDC_=false;
  }
  
  InitializeParameterFiles();
  InitializeHistograms();
  hddaq::unpacker::UnpackerManager& gUnpacker
    = hddaq::unpacker::GUnpacker::get_instance();

  gUnpacker.set_config_file( m_UnpackerFileName, m_DigitFileName, m_CmapFileName);

  return true;
}

bool ConfMan::InitializeEvDisp( void )
{
  static const std::string funcname = "[ConfMan::InitializeEvDisp]";
  evDisp_ = & EvDisp::GetInstance();
  evDisp_->Initialize();
}
