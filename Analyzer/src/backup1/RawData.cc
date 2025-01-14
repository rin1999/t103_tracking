/*
  RawData.cc

  2016/2  K.Shirotori
*/

#include "RawData.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "Random/Randomize.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cfloat>

#include <TROOT.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>

#include "DataType.hh"
#include "DetectorID.hh"
#include "PrimInfo.hh"
#include "HodoRawHit.hh"
#include "TrRawHit.hh"
#include "s_BeamRawHit.hh"
#include "s_ScatRawHit.hh"
#include "TemplateLib.hh"
#include "EventParameter.hh"

#include "ConfMan.hh"
#include "TrGeomMan.hh"

#define check1 0
#define check2 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
  T0RHC(0), RPCRHC(0),
  bSSDRHC(), sSSD1RHC(), sSSD2RHC(),
  s_BeamRHC(), s_ScatRHC()
{}

RawData::~RawData()
{
  clearAll();
}

bool RawData::AddTrRHit( TrRHitContainer& cont,
			 int Layer, int Wire, 
			 double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddTrRHit]";
 
  TrRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    TrRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
	q->WireId()==Wire ){
      p=q; break;
    }
  }
  if(!p){
    p = new TrRawHit( Layer, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetPosX( PosX );
    p->SetPosY( PosY );
    p->SetDL( DL );
    
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int Seg,
			   double TDC0_t, double TDC0_tot,
			   double TDC1_t, double TDC1_tot,
			   double ADC0_t, double ADC0_hgt,
			   double ADC1_t, double ADC1_hgt)
{
  static const std::string funcname = "[RawData::AddHodoRHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->SegmentId()==Seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetTdc0Time( TDC0_t );
    p->SetTdc0Tot( TDC0_tot );
    p->SetTdc1Time( TDC1_t );
    p->SetTdc1Tot( TDC1_tot );
    p->SetAdc0Time( ADC0_t );
    p->SetAdc0Hgt( ADC0_hgt );
    p->SetAdc1Time( ADC1_t );
    p->SetAdc1Hgt( ADC1_hgt );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

//For simple tracking
bool RawData::AddsBTrRHit( s_BeamRHitContainer& cont,
			   int TrackID, int Type,
			   int Layer, int Wire, 
			   double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddsBeamTrRHit]";
  
  s_BeamRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_BeamRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_BeamRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsTrRHit( Layer, Wire, PosX, PosY, DL );
  
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsBHodoRHit( s_BeamRHitContainer& cont,
			     int TrackID, int Type,
			     int DetId, int Seg,
			     double TDC0_t, double TDC0_tot,
			     double TDC1_t, double TDC1_tot,
			     double ADC0_t, double ADC0_hgt,
			     double ADC1_t, double ADC1_hgt )
{
  static const std::string funcname = "[RawData::AddsBHodoRHit]";

  s_BeamRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_BeamRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_BeamRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsHodoRHit( DetId, Seg,
		     TDC0_t, TDC0_tot, TDC1_t, TDC1_tot,
		     ADC0_t, ADC0_hgt, ADC1_t, ADC1_hgt );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsSTrRHit( s_ScatRHitContainer& cont,
			   int TrackID, int Type,
			   int Layer, int Wire, 
			   double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddsScatTrRHit]";
  
  s_ScatRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_ScatRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_ScatRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsTrRHit( Layer, Wire, PosX, PosY, DL );
  
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsSHodoRHit( s_ScatRHitContainer& cont,
			     int TrackID, int Type,
			     int DetId, int Seg,
			     double TDC0_t, double TDC0_tot,
			     double TDC1_t, double TDC1_tot,
			     double ADC0_t, double ADC0_hgt,
			     double ADC1_t, double ADC1_hgt )
{
  static const std::string funcname = "[RawData::AddsSHodoRHit]";

  s_ScatRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_ScatRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_ScatRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsHodoRHit( DetId, Seg,
		     TDC0_t, TDC0_tot, TDC1_t, TDC1_tot,
		     ADC0_t, ADC0_hgt, ADC1_t, ADC1_hgt );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

void RawData::clearAll()
{
  std::for_each(T0RHC.begin(), T0RHC.end(), DeleteObject());
  T0RHC.clear();

  std::for_each(RPCRHC.begin(), RPCRHC.end(), DeleteObject());
  RPCRHC.clear();

  for( int l=0; l<=PlMaxbSSD; ++l){
    for_each( bSSDRHC[l].begin(),  bSSDRHC[l].end(), DeleteObject());
    bSSDRHC[l].clear();
  }

  for( int l=0; l<=PlMaxsSSD1; ++l){
    for_each( sSSD1RHC[l].begin(),  sSSD1RHC[l].end(), DeleteObject());
    sSSD1RHC[l].clear();
  }

  for( int l=0; l<=PlMaxsSSD2; ++l){
    for_each( sSSD2RHC[l].begin(),  sSSD2RHC[l].end(), DeleteObject());
    sSSD2RHC[l].clear();
  }
    
  std::for_each(s_BeamRHC.begin(), s_BeamRHC.end(), DeleteObject());
  s_BeamRHC.clear();

  std::for_each(s_ScatRHC.begin(), s_ScatRHC.end(), DeleteObject());
  s_ScatRHC.clear();

  return;
}

bool RawData::DecodeRawHits( TFile* iFile, int i_entry )
{
  clearAll();
  ConfMan *confMan = ConfMan::GetConfManager();
  const TrGeomMan & geomMan=TrGeomMan::GetInstance();

  AddTrRHit(bSSDRHC[0], 0, 0, 0., 0., 0.);
  AddTrRHit(sSSD1RHC[0], 0, 0, 0., 0., 0.);
  AddTrRHit(sSSD2RHC[0], 0, 0, 0., 0., 0.);
 

  // Setup for input SSD data
  // const std::vector<std::vector<std::vector<int>>> SSD_pos ={
  //   {{0}, {1}},                    // SSD station 0
  //   {{2}, {3}},                    // SSD station 1
  //   {{4}, {5}, {6}},               // SSD station 2
  //   {{7}, {8}, {9}},               // SSD station 3
  //   {{10}, {11}},                  // SSD station 4
  //   {{12, 13}, {14, 15}, {16, 17}},// SSD station 5
  //   {{18, 19}, {20, 21}, {22, 23}},// SSD station 6
  //   {{24, 25}, {26, 27}}           // SSD station 7
  // };
  const std::vector<std::vector<int>> SSD_pos ={
    {0, 1},      //SSD station 0
    {2, 3},      //SSD station 1
    {4, 5, 6},   //SSD station 2
    {7, 8, 9},   //SSD station 3
    {10, 11},    //SSD station 4
    {12, 13, 14},//SSD station 5
    {15, 16, 17},//SSD station 6
    {18, 19}     //SSD station 7
  };

  TDirectoryFile *iDir = (TDirectoryFile*)iFile->Get("out1");
  TTree *tree = (TTree*)iDir->Get("T0AnaTree");

  std::vector<std::vector<double>> *tdc_top_lead = 0;
  std::vector<std::vector<double>> *tdc_top_tot = 0;
  std::vector<std::vector<double>> *tdc_bot_lead = 0;
  std::vector<std::vector<double>> *tdc_bot_tot = 0;

  std::array<double, NumOfSegT0> adc_top_hgt;
  std::array<double, NumOfSegT0> adc_top_t;
  std::array<double, NumOfSegT0> adc_bot_hgt;
  std::array<double, NumOfSegT0> adc_bot_t;

  std::vector<std::vector<double>> *rpc_lft_lead = 0;
  std::vector<std::vector<double>> *rpc_lft_tot = 0;
  std::vector<std::vector<double>> *rpc_rgt_lead = 0;
  std::vector<std::vector<double>> *rpc_rgt_tot = 0;


  std::vector<double> *bSSD_cl_station = 0;
  std::vector<double> *bSSD_cl_plane = 0;
  std::vector<double> *bSSD_cl_sens = 0;
  std::vector<double> *bSSD_cl_wgtavgstrip = 0;

  std::vector<double> *bSSD_ls_x0_x = 0;
  std::vector<double> *bSSD_ls_x0_y = 0;
  std::vector<double> *bSSD_ls_x0_z = 0;
  std::vector<double> *bSSD_ls_x1_x = 0;
  std::vector<double> *bSSD_ls_x1_y = 0;
  std::vector<double> *bSSD_ls_x1_z = 0;
  

  std::vector<double> *sSSD1_cl_station = 0;
  std::vector<double> *sSSD1_cl_plane = 0;
  std::vector<double> *sSSD1_cl_sens = 0;
  std::vector<double> *sSSD1_cl_wgtavgstrip = 0;

  std::vector<double> *sSSD1_ls_x0_x = 0;
  std::vector<double> *sSSD1_ls_x0_y = 0;
  std::vector<double> *sSSD1_ls_x0_z = 0;
  std::vector<double> *sSSD1_ls_x1_x = 0;
  std::vector<double> *sSSD1_ls_x1_y = 0;
  std::vector<double> *sSSD1_ls_x1_z = 0;
  

  std::vector<double> *sSSD2_cl_station = 0;
  std::vector<double> *sSSD2_cl_plane = 0;
  std::vector<double> *sSSD2_cl_sens = 0;
  std::vector<double> *sSSD2_cl_wgtavgstrip = 0;

  std::vector<double> *sSSD2_ls_x0_x = 0;
  std::vector<double> *sSSD2_ls_x0_y = 0;
  std::vector<double> *sSSD2_ls_x0_z = 0;
  std::vector<double> *sSSD2_ls_x1_x = 0;
  std::vector<double> *sSSD2_ls_x1_y = 0;
  std::vector<double> *sSSD2_ls_x1_z = 0;
  

  tree->SetBranchAddress("TDC_top_lead",  &tdc_top_lead);
  tree->SetBranchAddress("TDC_top_tot",  &tdc_top_tot);
  tree->SetBranchAddress("TDC_bot_lead",  &tdc_bot_lead);
  tree->SetBranchAddress("TDC_bot_tot",  &tdc_bot_tot);

  tree->SetBranchAddress("ADC_top_hgt", &adc_top_hgt);
  tree->SetBranchAddress("ADC_top_t", &adc_top_t);
  tree->SetBranchAddress("ADC_bot_hgt", &adc_bot_hgt);
  tree->SetBranchAddress("ADC_bot_t", &adc_bot_t);

  tree->SetBranchAddress("RPC_lft_lead",  &rpc_lft_lead);
  tree->SetBranchAddress("RPC_lft_tot",  &rpc_lft_tot);
  tree->SetBranchAddress("RPC_rgt_lead",  &rpc_rgt_lead);
  tree->SetBranchAddress("RPC_rgt_tot",  &rpc_rgt_tot);


  tree->SetBranchAddress("bSSD_cl_station", &bSSD_cl_station);
  tree->SetBranchAddress("bSSD_cl_plane",   &bSSD_cl_plane);
  tree->SetBranchAddress("bSSD_cl_sens",    &bSSD_cl_sens);
  tree->SetBranchAddress("bSSD_cl_wgtavgstrip", &bSSD_cl_wgtavgstrip);

  tree->SetBranchAddress("bSSD_ls_x0_x", &bSSD_ls_x0_x);
  tree->SetBranchAddress("bSSD_ls_x0_y", &bSSD_ls_x0_y);
  tree->SetBranchAddress("bSSD_ls_x0_z", &bSSD_ls_x0_z);
  tree->SetBranchAddress("bSSD_ls_x1_x", &bSSD_ls_x1_x);
  tree->SetBranchAddress("bSSD_ls_x1_y", &bSSD_ls_x1_y);
  tree->SetBranchAddress("bSSD_ls_x1_z", &bSSD_ls_x1_z);


  tree->SetBranchAddress("sSSD1_cl_station", &sSSD1_cl_station);
  tree->SetBranchAddress("sSSD1_cl_plane",   &sSSD1_cl_plane);
  tree->SetBranchAddress("sSSD1_cl_sens",    &sSSD1_cl_sens);
  tree->SetBranchAddress("sSSD1_cl_wgtavgstrip", &sSSD1_cl_wgtavgstrip);

  tree->SetBranchAddress("sSSD1_ls_x0_x", &sSSD1_ls_x0_x);
  tree->SetBranchAddress("sSSD1_ls_x0_y", &sSSD1_ls_x0_y);
  tree->SetBranchAddress("sSSD1_ls_x0_z", &sSSD1_ls_x0_z);
  tree->SetBranchAddress("sSSD1_ls_x1_x", &sSSD1_ls_x1_x);
  tree->SetBranchAddress("sSSD1_ls_x1_y", &sSSD1_ls_x1_y);
  tree->SetBranchAddress("sSSD1_ls_x1_z", &sSSD1_ls_x1_z);


  tree->SetBranchAddress("sSSD2_cl_station", &sSSD2_cl_station);
  tree->SetBranchAddress("sSSD2_cl_plane",   &sSSD2_cl_plane);
  tree->SetBranchAddress("sSSD2_cl_sens",    &sSSD2_cl_sens);
  tree->SetBranchAddress("sSSD2_cl_wgtavgstrip", &sSSD2_cl_wgtavgstrip);

  tree->SetBranchAddress("sSSD2_ls_x0_x", &sSSD2_ls_x0_x);
  tree->SetBranchAddress("sSSD2_ls_x0_y", &sSSD2_ls_x0_y);
  tree->SetBranchAddress("sSSD2_ls_x0_z", &sSSD2_ls_x0_z);
  tree->SetBranchAddress("sSSD2_ls_x1_x", &sSSD2_ls_x1_x);
  tree->SetBranchAddress("sSSD2_ls_x1_y", &sSSD2_ls_x1_y);
  tree->SetBranchAddress("sSSD2_ls_x1_z", &sSSD2_ls_x1_z);


  tree->GetEntry(i_entry);


  // Loop for T0 signals
  for(int i_seg_t0 = 0; i_seg_t0 < NumOfSegT0; i_seg_t0++){
    double t0_tdc_top_t = -DBL_MAX;
    double t0_tdc_top_tot = -DBL_MAX;
    double t0_adc_top_t = -DBL_MAX;
    double t0_adc_top_hgt = -DBL_MAX;

    int n_vec_t0_top = tdc_top_lead->at(i_seg_t0).size();
    for(int i_vec_t0_top = 0; i_vec_t0_top < n_vec_t0_top; i_vec_t0_top++){
      if(tdc_top_lead->at(i_seg_t0).at(i_vec_t0_top) >= T0TdcMin && tdc_top_lead->at(i_seg_t0).at(i_vec_t0_top) <= T0TdcMax && !tdc_top_tot->at(i_seg_t0).empty()){
	t0_tdc_top_t = tdc_top_lead->at(i_seg_t0).at(i_vec_t0_top);
	t0_tdc_top_tot = tdc_top_tot->at(i_seg_t0).at(0);

	t0_adc_top_t = adc_top_t.at(i_seg_t0);
	t0_adc_top_hgt = adc_top_hgt.at(i_seg_t0)*ChToADC;

	break;// Get out of the for loop at fitst hit
      }//if(tdc_gate)
    }//for(i_vec_t0_top:n_vec_t0_top)


    double t0_tdc_bot_t = -DBL_MAX;
    double t0_tdc_bot_tot = -DBL_MAX;
    double t0_adc_bot_t = -DBL_MAX;
    double t0_adc_bot_hgt = -DBL_MAX;

    int n_vec_t0_bot = tdc_bot_lead->at(i_seg_t0).size();
    for(int i_vec_t0_bot = 0; i_vec_t0_bot < n_vec_t0_bot; i_vec_t0_bot++){
      if(tdc_bot_lead->at(i_seg_t0).at(i_vec_t0_bot) >= T0TdcMin && tdc_bot_lead->at(i_seg_t0).at(i_vec_t0_bot) <= T0TdcMax && !tdc_bot_tot->at(i_seg_t0).empty()){
	t0_tdc_bot_t = tdc_bot_lead->at(i_seg_t0).at(i_vec_t0_bot);
	t0_tdc_bot_tot = tdc_bot_tot->at(i_seg_t0).at(0);

	t0_adc_bot_t = adc_bot_t.at(i_seg_t0);
	t0_adc_bot_hgt = adc_bot_hgt.at(i_seg_t0)*ChToADC;

	break;// Get out of the for loop at fitst hit
      }//if(tdc_gate)
    }//for(i_vec_t0_bot:n_vec_t0_bot)

    if(t0_tdc_top_t == -DBL_MAX || t0_tdc_bot_t == -DBL_MAX){
      t0_tdc_top_t = -DBL_MAX;
      t0_tdc_top_tot = -DBL_MAX;
      t0_adc_top_t = -DBL_MAX;
      t0_adc_top_hgt = -DBL_MAX;

      t0_tdc_bot_t = -DBL_MAX;
      t0_tdc_bot_tot = -DBL_MAX;
      t0_adc_bot_t = -DBL_MAX;
      t0_adc_bot_hgt = -DBL_MAX;
    }//if(Reset parameter at no hit)

    AddHodoRHit(T0RHC,
		IdOfT0, i_seg_t0,
		t0_tdc_top_t, t0_tdc_top_tot,
		t0_tdc_bot_t, t0_tdc_bot_tot,
		t0_adc_top_t, t0_adc_top_hgt,
		t0_adc_bot_t, t0_adc_bot_hgt);
    AddsBHodoRHit(s_BeamRHC, 0, 0,
		  IdOfT0, i_seg_t0,
		  t0_tdc_top_t, t0_tdc_top_tot,
		  t0_tdc_bot_t, t0_tdc_bot_tot,
		  t0_adc_top_t, t0_adc_top_hgt,
		  t0_adc_bot_t, t0_adc_bot_hgt);
  }//for(i_seg_t0:NumOfSegT0)


  // Loop for RPC signals
  for(int i_seg_rpc = 0; i_seg_rpc < NumOfSegRPC; i_seg_rpc++){
    double rpc_tdc_lft_t = -DBL_MAX;
    double rpc_tdc_lft_tot = -DBL_MAX;
    double rpc_adc_lft_t = -DBL_MAX;
    double rpc_adc_lft_hgt = -DBL_MAX;

    int n_vec_rpc_lft = rpc_lft_lead->at(i_seg_rpc).size();
    for(int i_vec_rpc_lft = 0; i_vec_rpc_lft < n_vec_rpc_lft; i_vec_rpc_lft++){
      if(rpc_lft_lead->at(i_seg_rpc).at(i_vec_rpc_lft) >= RPCTdcMin && rpc_lft_lead->at(i_seg_rpc).at(i_vec_rpc_lft) <= RPCTdcMax && !rpc_lft_tot->at(i_seg_rpc).empty()){
	rpc_tdc_lft_t = rpc_lft_lead->at(i_seg_rpc).at(i_vec_rpc_lft);
	rpc_tdc_lft_tot = rpc_lft_tot->at(i_seg_rpc).at(0);

	break;// Get out of the for loop at fitst hit
      }//if(tdc_gate)
    }//for(i_vec_rpc_lft:n_vec_rpc_lft)


    double rpc_tdc_rgt_t = -DBL_MAX;
    double rpc_tdc_rgt_tot = -DBL_MAX;
    double rpc_adc_rgt_t = -DBL_MAX;
    double rpc_adc_rgt_hgt = -DBL_MAX;

    int n_vec_rpc_rgt = rpc_rgt_lead->at(i_seg_rpc).size();
    for(int i_vec_rpc_rgt = 0; i_vec_rpc_rgt < n_vec_rpc_rgt; i_vec_rpc_rgt++){
      if(rpc_rgt_lead->at(i_seg_rpc).at(i_vec_rpc_rgt) >= RPCTdcMin && rpc_rgt_lead->at(i_seg_rpc).at(i_vec_rpc_rgt) <= RPCTdcMax && !rpc_rgt_tot->at(i_seg_rpc).empty()){
	rpc_tdc_rgt_t = rpc_rgt_lead->at(i_seg_rpc).at(i_vec_rpc_rgt);
	rpc_tdc_rgt_tot = rpc_rgt_tot->at(i_seg_rpc).at(0);

	break;// Get out of the for loop at fitst hit
      }//if(tdc_gate)
    }//for(i_vec_rpc_rgt:n_vec_rpc_rgt)

    if(rpc_tdc_lft_t == -DBL_MAX || rpc_tdc_rgt_t == -DBL_MAX){
      rpc_tdc_lft_t = -DBL_MAX;
      rpc_tdc_lft_tot = -DBL_MAX;
      rpc_adc_lft_t = -DBL_MAX;
      rpc_adc_lft_hgt = -DBL_MAX;

      rpc_tdc_rgt_t = -DBL_MAX;
      rpc_tdc_rgt_tot = -DBL_MAX;
      rpc_adc_rgt_t = -DBL_MAX;
      rpc_adc_rgt_hgt = -DBL_MAX;
    }

    AddHodoRHit(RPCRHC,
		IdOfRPC, i_seg_rpc,
		rpc_tdc_lft_t, rpc_tdc_lft_tot,
		rpc_tdc_rgt_t, rpc_tdc_rgt_tot,
		rpc_adc_lft_t, rpc_adc_lft_hgt,
		rpc_adc_rgt_t, rpc_adc_rgt_hgt);
    AddsSHodoRHit(s_ScatRHC, 0, 0,
		  IdOfRPC, i_seg_rpc,
		  rpc_tdc_lft_t, rpc_tdc_lft_tot,
		  rpc_tdc_rgt_t, rpc_tdc_rgt_tot,
		  rpc_adc_lft_t, rpc_adc_lft_hgt,
		  rpc_adc_rgt_t, rpc_adc_rgt_hgt);
  }//for(i_seg_rpc:NumOfSegRPC)


  // Loop for bSSD signals
  int n_vec_bssd = bSSD_cl_station->size();
  for(int i_vec = 0; i_vec < n_vec_bssd; i_vec++){
    int bSSD_cl_pos = SSD_pos[bSSD_cl_station->at(i_vec)][bSSD_cl_plane->at(i_vec)] + 1;
    double bSSD_ls_x = (bSSD_ls_x0_x->at(i_vec) + bSSD_ls_x1_x->at(i_vec))/2.0;
    double bSSD_ls_y = (bSSD_ls_x0_y->at(i_vec) + bSSD_ls_x1_y->at(i_vec))/2.0;
    double angle = geomMan.GetTiltAngle( bSSD_cl_pos );
    double l = bSSD_ls_x*cos(angle*Deg2Rad) + bSSD_ls_y*sin(angle*Deg2Rad);
    double dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetbSSDResol() );

    AddTrRHit(bSSDRHC[bSSD_cl_pos-PlOffsbSSD], bSSD_cl_pos,
	      0, bSSD_ls_x, bSSD_ls_y, dl);
    AddsBTrRHit(s_BeamRHC, 0, 0, bSSD_cl_pos,
		0, bSSD_ls_x, bSSD_ls_y, dl);
  }//for(i_vec:n_vec_bssd)


  // Loop for sSSD1 signals
  int n_vec_sssd1 = sSSD1_cl_station->size();
  for(int i_vec = 0; i_vec < n_vec_sssd1; i_vec++){
    int sSSD1_cl_pos = SSD_pos[sSSD1_cl_station->at(i_vec)][sSSD1_cl_plane->at(i_vec)] + 1;
    double sSSD1_ls_x = (sSSD1_ls_x0_x->at(i_vec) + sSSD1_ls_x1_x->at(i_vec))/2.0;
    double sSSD1_ls_y = (sSSD1_ls_x0_y->at(i_vec) + sSSD1_ls_x1_y->at(i_vec))/2.0;
    double angle = geomMan.GetTiltAngle( sSSD1_cl_pos );
    double l = sSSD1_ls_x*cos(angle*Deg2Rad) + sSSD1_ls_y*sin(angle*Deg2Rad);
    double dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetsSSD1Resol() );

    AddTrRHit(sSSD1RHC[sSSD1_cl_pos-PlOffssSSD1], sSSD1_cl_pos,
	      0, sSSD1_ls_x, sSSD1_ls_y, dl);
    AddsSTrRHit(s_ScatRHC, 0, 0, sSSD1_cl_pos,
		0, sSSD1_ls_x, sSSD1_ls_y, dl);
  }//for(i_vec:n_vec_sssd2)


  // Loop for sSSD2 signals
  int n_vec_sssd2 = sSSD2_cl_station->size();
  for(int i_vec = 0; i_vec < n_vec_sssd2; i_vec++){
    int sSSD2_cl_pos = SSD_pos[sSSD2_cl_station->at(i_vec)][sSSD2_cl_plane->at(i_vec)] + 1;
    double sSSD2_ls_x = (sSSD2_ls_x0_x->at(i_vec) + sSSD2_ls_x1_x->at(i_vec))/2.0;
    double sSSD2_ls_y = (sSSD2_ls_x0_y->at(i_vec) + sSSD2_ls_x1_y->at(i_vec))/2.0;
    double angle = geomMan.GetTiltAngle( sSSD2_cl_pos );
    double l = sSSD2_ls_x*cos(angle*Deg2Rad) + sSSD2_ls_y*sin(angle*Deg2Rad);
    double dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetsSSD2Resol() );

    if(sSSD2_cl_sens->at(i_vec)==0){
      if(sSSD2_ls_x <= 0.0 && sSSD2_ls_y <= 0.0){
	AddTrRHit(sSSD2RHC[sSSD2_cl_pos-PlOffssSSD2], sSSD2_cl_pos,
		  0, sSSD2_ls_x, sSSD2_ls_y, dl);
	AddsSTrRHit(s_ScatRHC, 0, 0, sSSD2_cl_pos,
		    0, sSSD2_ls_x, sSSD2_ls_y, dl);
      }//Position selection
    }else{
      if(sSSD2_ls_x >= 0.0 && sSSD2_ls_y >= 0.0){
	AddTrRHit(sSSD2RHC[sSSD2_cl_pos-PlOffssSSD2], sSSD2_cl_pos,
		  0, sSSD2_ls_x, sSSD2_ls_y, dl);
	AddsSTrRHit(s_ScatRHC, 0, 0, sSSD2_cl_pos,
		    0, sSSD2_ls_x, sSSD2_ls_y, dl);
      }//Position selection
    }//if(sSSD2_cl_sens)
  }//for(i_vec:n_vec_sssd1)

  return true;
}// RawData::DecodeRawHits

const HodoRHitContainer& RawData::GetT0RHC() const
{
  return T0RHC;
}

const HodoRHitContainer& RawData::GetRPCRHC() const
{
  return RPCRHC;
}

const TrRHitContainer & RawData::GetbSSDRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxbSSD ) layer=0;
  return bSSDRHC[layer];
}

const TrRHitContainer & RawData::GetsSSD1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxsSSD1 ) layer=0;
  return sSSD1RHC[layer];
}

const TrRHitContainer & RawData::GetsSSD2RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxsSSD2 ) layer=0;
  return sSSD2RHC[layer];
}

const s_BeamRHitContainer & RawData::GetsBeamRHC( void ) const
{
  return s_BeamRHC;
}

const s_ScatRHitContainer & RawData::GetsScatRHC( void ) const
{
  return s_ScatRHC;
}
