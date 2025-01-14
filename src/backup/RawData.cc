/*
  RawData.cc

  2024/04  K.Shirotori
*/

#include "RawData.hh"

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

#include "DetectorInfo.hh"
#include "DCRawHit.hh"
#include "TrRawHit.hh"
#include "HodoRawHit.hh"

#include "TemplateLib.hh"

#include "ConfMan.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
   BDCRHC(),
   KLDCRHC(),
   BFTRHC(),
   SFTRHC(),

   UTOFRHC(0),
   DTOFRHC(0),
   LTOFRHC(0),
   T0RHC(0),
   T0rRHC(0),
   BrefRHC(0),
   T1RHC(0),
   BHTRHC(0)
{}

RawData::~RawData()
{
  clearAll();
}

bool RawData::AddDCRHit( DCRHitContainer& cont,
			 int LayerId, int WireId, 
			 DoubleVec lTdc, DoubleVec Tot )
{
   static const std::string funcname = "[RawData::AddDCRHit]";
   
   DCRawHit *p=0;
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      DCRawHit *q=cont[i];
      if( q->LayerId()==LayerId &&
          q->WireId()==WireId ){
         p=q; break;
      }
   }
   if(!p){
      p = new DCRawHit( LayerId, WireId );
      if(p) cont.push_back(p);
   }
   if(p){
      p->SetlTdc( lTdc );
      p->SetTot( Tot );
      
      return true;
   }else{
      std::cerr << funcname << ": new fail." << std::endl;
      return false;
   }
}

bool RawData::AddTrRHit( TrRHitContainer& cont,
			 int LayerId, int FiberId, 
			 DoubleVec lTdc, DoubleVec Tot )
{
   static const std::string funcname = "[RawData::AddDCRHit]";
   
   TrRawHit *p=0;
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      TrRawHit *q=cont[i];
      if( q->LayerId()==LayerId &&
          q->FiberId()==FiberId ){
         p=q; break;
      }
   }
   if(!p){
      p = new TrRawHit( LayerId, FiberId );
      if(p) cont.push_back(p);
   }
   if(p){
      p->SetlTdc( lTdc );
      p->SetTot( Tot );
      
      return true;
   }else{
      std::cerr << funcname << ": new fail." << std::endl;
      return false;
   }
}

bool RawData::AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int LayId, int SegId,
			   DoubleVec lTdc1, DoubleVec Tot1,
			   DoubleVec lTdc2, DoubleVec Tot2)
{
   static const std::string funcname = "[RawData::AddHodoRHit]";
   
   HodoRawHit *p=0;
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *q=cont[i];
      if( q->DetectorId()==DetId &&
          q->LayerId()==LayId &&
          q->SegmentId()==SegId){
         p=q; break;
      }
   }
   if(!p){
      p = new HodoRawHit( DetId, LayId, SegId );
      if(p) cont.push_back(p);
   }
   if(p){
      p->SetlTdc1( lTdc1 );
      p->SetlTdc2( lTdc2 );
      p->SetTot1( Tot1 );
      p->SetTot2( Tot2 );

      return true;
   }else{
      std::cerr << funcname << ": new fail." << std::endl;
      return false;
   }
}

void RawData::clearAll()
{
   for( int l=0; l<=NumOfLayersBDC; ++l){
      std::for_each( BDCRHC[l].begin(),  BDCRHC[l].end(), DeleteObject());
      BDCRHC[l].clear();
   }

   for( int l=0; l<=NumOfLayersKLDC; ++l){
      std::for_each( KLDCRHC[l].begin(),  KLDCRHC[l].end(), DeleteObject());
      KLDCRHC[l].clear();
   }

   for( int l=0; l<=NumOfLayersBFT; ++l){
      std::for_each( BFTRHC[l].begin(),  BFTRHC[l].end(), DeleteObject());
      BFTRHC[l].clear();
   }

   for( int l=0; l<=NumOfLayersSFT; ++l){
      std::for_each( SFTRHC[l].begin(),  SFTRHC[l].end(), DeleteObject());
      SFTRHC[l].clear();
   }

   std::for_each(UTOFRHC.begin(), UTOFRHC.end(), DeleteObject());
   UTOFRHC.clear();

   std::for_each(DTOFRHC.begin(), DTOFRHC.end(), DeleteObject());
   DTOFRHC.clear();
   
   std::for_each(LTOFRHC.begin(), LTOFRHC.end(), DeleteObject());
   LTOFRHC.clear();

   std::for_each(T0RHC.begin(), T0RHC.end(), DeleteObject());
   T0RHC.clear();

   std::for_each(T0rRHC.begin(), T0rRHC.end(), DeleteObject());
   T0rRHC.clear();
   
   std::for_each(BrefRHC.begin(), BrefRHC.end(), DeleteObject());
   BrefRHC.clear();

   std::for_each(T1RHC.begin(), T1RHC.end(), DeleteObject());
   T1RHC.clear();

   std::for_each(BHTRHC.begin(), BHTRHC.end(), DeleteObject());
   BHTRHC.clear();

   return;
}

bool RawData::DecodeRawHits( TFile* iFile, int i_entry )
{
  clearAll();

  TTree *tree = (TTree*)iFile->Get("tree");
  
  std::vector<std::vector<double>> *tot_utof_l = 0;
  std::vector<std::vector<double>> *tot_utof_r = 0;
  std::vector<std::vector<double>> *tot_dtof_l = 0;
  std::vector<std::vector<double>> *tot_dtof_r = 0;
  std::vector<std::vector<double>> *tot_ltof_l = 0;
  std::vector<std::vector<double>> *tot_ltof_r = 0;
  std::vector<std::vector<double>> *tot_t0_l = 0;
  std::vector<std::vector<double>> *tot_t0_r = 0;
  std::vector<std::vector<double>> *tot_t0r_l = 0;
  std::vector<std::vector<double>> *tot_t0r_r = 0;
  std::vector<std::vector<double>> *tot_bref_l = 0;
  std::vector<std::vector<double>> *tot_bref_r = 0;
  std::vector<std::vector<double>> *tot_t1_l = 0;
  std::vector<std::vector<double>> *tot_t1_r = 0;
  std::vector<std::vector<double>> *tot_bht_l = 0;
  std::vector<std::vector<double>> *tot_bht_r = 0;

  std::vector<std::vector<double>> *ltdc_utof_l = 0;
  std::vector<std::vector<double>> *ltdc_utof_r = 0;
  std::vector<std::vector<double>> *ltdc_dtof_l = 0;
  std::vector<std::vector<double>> *ltdc_dtof_r = 0;
  std::vector<std::vector<double>> *ltdc_ltof_l = 0;
  std::vector<std::vector<double>> *ltdc_ltof_r = 0;
  std::vector<std::vector<double>> *ltdc_t0_l = 0;
  std::vector<std::vector<double>> *ltdc_t0_r = 0;
  std::vector<std::vector<double>> *ltdc_t0r_l = 0;
  std::vector<std::vector<double>> *ltdc_t0r_r = 0;
  std::vector<std::vector<double>> *ltdc_bref_l = 0;
  std::vector<std::vector<double>> *ltdc_bref_r = 0;
  std::vector<std::vector<double>> *ltdc_t1_l = 0;
  std::vector<std::vector<double>> *ltdc_t1_r = 0;
  std::vector<std::vector<double>> *ltdc_bht_l = 0;
  std::vector<std::vector<double>> *ltdc_bht_r = 0;

  std::vector<std::vector<double>> *tot_bdc_l1 = 0;
  std::vector<std::vector<double>> *tot_bdc_l2 = 0;
  std::vector<std::vector<double>> *tot_bdc_l3 = 0;
  std::vector<std::vector<double>> *tot_bdc_l4 = 0;
  std::vector<std::vector<double>> *tot_bdc_l5 = 0;
  std::vector<std::vector<double>> *tot_bdc_l6 = 0;
  std::vector<std::vector<double>> *tot_bdc_l7 = 0;
  std::vector<std::vector<double>> *tot_bdc_l8 = 0;

  std::vector<std::vector<double>> *ltdc_bdc_l1 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l2 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l3 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l4 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l5 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l6 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l7 = 0;
  std::vector<std::vector<double>> *ltdc_bdc_l8 = 0;

  std::vector<std::vector<double>> *tot_kldc_l1 = 0;
  std::vector<std::vector<double>> *tot_kldc_l2 = 0;
  std::vector<std::vector<double>> *tot_kldc_l3 = 0;
  std::vector<std::vector<double>> *tot_kldc_l4 = 0;
  std::vector<std::vector<double>> *tot_kldc_l5 = 0;
  std::vector<std::vector<double>> *tot_kldc_l6 = 0;
  std::vector<std::vector<double>> *tot_kldc_l7 = 0;
  std::vector<std::vector<double>> *tot_kldc_l8 = 0;

  std::vector<std::vector<double>> *ltdc_kldc_l1 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l2 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l3 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l4 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l5 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l6 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l7 = 0;
  std::vector<std::vector<double>> *ltdc_kldc_l8 = 0;

  std::vector<std::vector<double>> *tot_bft_l1 = 0;
  std::vector<std::vector<double>> *tot_bft_l2 = 0;
  std::vector<std::vector<double>> *tot_bft_l3 = 0;
  std::vector<std::vector<double>> *tot_bft_l4 = 0;
  std::vector<std::vector<double>> *tot_bft_l5 = 0;
  std::vector<std::vector<double>> *tot_bft_l6 = 0;

  std::vector<std::vector<double>> *ltdc_bft_l1 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l2 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l3 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l4 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l5 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l6 = 0;

  std::vector<std::vector<double>> *tot_sft_l1 = 0;
  std::vector<std::vector<double>> *tot_sft_l2 = 0;
  std::vector<std::vector<double>> *tot_sft_l3 = 0;
  std::vector<std::vector<double>> *tot_sft_l4 = 0;
  std::vector<std::vector<double>> *tot_sft_l5 = 0;
  std::vector<std::vector<double>> *tot_sft_l6 = 0;

  std::vector<std::vector<double>> *ltdc_sft_l1 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l2 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l3 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l4 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l5 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l6 = 0;

  TBranch *Btot_utof_l = 0;
  TBranch *Btot_utof_r = 0;
  TBranch *Btot_dtof_l = 0;
  TBranch *Btot_dtof_r = 0;
  TBranch *Btot_ltof_l = 0;
  TBranch *Btot_ltof_r = 0;
  TBranch *Btot_t0_l = 0;
  TBranch *Btot_t0_r = 0;
  TBranch *Btot_t0r_l = 0;
  TBranch *Btot_t0r_r = 0;
  TBranch *Btot_bref_l = 0;
  TBranch *Btot_bref_r = 0;
  TBranch *Btot_t1_l = 0;
  TBranch *Btot_t1_r = 0;
  TBranch *Btot_bht_l = 0;
  TBranch *Btot_bht_r = 0;

  TBranch *Bltdc_utof_l = 0;
  TBranch *Bltdc_utof_r = 0;
  TBranch *Bltdc_dtof_l = 0;
  TBranch *Bltdc_dtof_r = 0;
  TBranch *Bltdc_ltof_l = 0;
  TBranch *Bltdc_ltof_r = 0;
  TBranch *Bltdc_t0_l = 0;
  TBranch *Bltdc_t0_r = 0;
  TBranch *Bltdc_t0r_l = 0;
  TBranch *Bltdc_t0r_r = 0;
  TBranch *Bltdc_bref_l = 0;
  TBranch *Bltdc_bref_r = 0;
  TBranch *Bltdc_t1_l = 0;
  TBranch *Bltdc_t1_r = 0;
  TBranch *Bltdc_bht_l = 0;
  TBranch *Bltdc_bht_r = 0;

  TBranch *Btot_bdc_l1 = 0;
  TBranch *Btot_bdc_l2 = 0;
  TBranch *Btot_bdc_l3 = 0;
  TBranch *Btot_bdc_l4 = 0;
  TBranch *Btot_bdc_l5 = 0;
  TBranch *Btot_bdc_l6 = 0;
  TBranch *Btot_bdc_l7 = 0;
  TBranch *Btot_bdc_l8 = 0;

  TBranch *Bltdc_bdc_l1 = 0;
  TBranch *Bltdc_bdc_l2 = 0;
  TBranch *Bltdc_bdc_l3 = 0;
  TBranch *Bltdc_bdc_l4 = 0;
  TBranch *Bltdc_bdc_l5 = 0;
  TBranch *Bltdc_bdc_l6 = 0;
  TBranch *Bltdc_bdc_l7 = 0;
  TBranch *Bltdc_bdc_l8 = 0;

  TBranch *Btot_kldc_l1 = 0;
  TBranch *Btot_kldc_l2 = 0;
  TBranch *Btot_kldc_l3 = 0;
  TBranch *Btot_kldc_l4 = 0;
  TBranch *Btot_kldc_l5 = 0;
  TBranch *Btot_kldc_l6 = 0;
  TBranch *Btot_kldc_l7 = 0;
  TBranch *Btot_kldc_l8 = 0;

  TBranch *Bltdc_kldc_l1 = 0;
  TBranch *Bltdc_kldc_l2 = 0;
  TBranch *Bltdc_kldc_l3 = 0;
  TBranch *Bltdc_kldc_l4 = 0;
  TBranch *Bltdc_kldc_l5 = 0;
  TBranch *Bltdc_kldc_l6 = 0;
  TBranch *Bltdc_kldc_l7 = 0;
  TBranch *Bltdc_kldc_l8 = 0;

  TBranch *Btot_bft_l1 = 0;
  TBranch *Btot_bft_l2 = 0;
  TBranch *Btot_bft_l3 = 0;
  TBranch *Btot_bft_l4 = 0;
  TBranch *Btot_bft_l5 = 0;
  TBranch *Btot_bft_l6 = 0;

  TBranch *Bltdc_bft_l1 = 0;
  TBranch *Bltdc_bft_l2 = 0;
  TBranch *Bltdc_bft_l3 = 0;
  TBranch *Bltdc_bft_l4 = 0;
  TBranch *Bltdc_bft_l5 = 0;
  TBranch *Bltdc_bft_l6 = 0;

  TBranch *Btot_sft_l1 = 0;
  TBranch *Btot_sft_l2 = 0;
  TBranch *Btot_sft_l3 = 0;
  TBranch *Btot_sft_l4 = 0;
  TBranch *Btot_sft_l5 = 0;
  TBranch *Btot_sft_l6 = 0;

  TBranch *Bltdc_sft_l1 = 0;
  TBranch *Bltdc_sft_l2 = 0;
  TBranch *Bltdc_sft_l3 = 0;
  TBranch *Bltdc_sft_l4 = 0;
  TBranch *Bltdc_sft_l5 = 0;
  TBranch *Bltdc_sft_l6 = 0;
  
  tree->SetBranchAddress("tot_utof_l", &tot_utof_l, &Btot_utof_l);
  tree->SetBranchAddress("tot_utof_r", &tot_utof_r, &Btot_utof_r);
  tree->SetBranchAddress("tot_dtof_l", &tot_dtof_l, &Btot_dtof_l);
  tree->SetBranchAddress("tot_dtof_r", &tot_dtof_r, &Btot_dtof_r);
  tree->SetBranchAddress("tot_ltof_l", &tot_ltof_l, &Btot_ltof_l);
  tree->SetBranchAddress("tot_ltof_r", &tot_ltof_r, &Btot_ltof_r);
  tree->SetBranchAddress("tot_t0_l", &tot_t0_l, &Btot_t0_l);
  tree->SetBranchAddress("tot_t0_r", &tot_t0_r, &Btot_t0_r);
  tree->SetBranchAddress("tot_t0r_l", &tot_t0r_l, &Btot_t0r_l);
  tree->SetBranchAddress("tot_t0r_r", &tot_t0r_r, &Btot_t0r_r);
  tree->SetBranchAddress("tot_bref_l", &tot_bref_l, &Btot_bref_l);
  tree->SetBranchAddress("tot_bref_r", &tot_bref_r, &Btot_bref_r);
  tree->SetBranchAddress("tot_t1_l", &tot_t1_l, &Btot_t1_l);
  tree->SetBranchAddress("tot_t1_r", &tot_t1_r, &Btot_t1_r);
  tree->SetBranchAddress("tot_bht_l", &tot_bht_l, &Btot_bht_l);
  tree->SetBranchAddress("tot_bht_r", &tot_bht_r, &Btot_bht_r);

  tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l);
  tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r);
  tree->SetBranchAddress("ltdc_dtof_l", &ltdc_dtof_l, &Bltdc_dtof_l);
  tree->SetBranchAddress("ltdc_dtof_r", &ltdc_dtof_r, &Bltdc_dtof_r);
  tree->SetBranchAddress("ltdc_ltof_l", &ltdc_ltof_l, &Bltdc_ltof_l);
  tree->SetBranchAddress("ltdc_ltof_r", &ltdc_ltof_r, &Bltdc_ltof_r);
  tree->SetBranchAddress("ltdc_t0_l", &ltdc_t0_l, &Bltdc_t0_l);
  tree->SetBranchAddress("ltdc_t0_r", &ltdc_t0_r, &Bltdc_t0_r);
  tree->SetBranchAddress("ltdc_t0r_l", &ltdc_t0r_l, &Bltdc_t0r_l);
  tree->SetBranchAddress("ltdc_t0r_r", &ltdc_t0r_r, &Bltdc_t0r_r);
  tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l);
  tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r);
  tree->SetBranchAddress("ltdc_t1_l", &ltdc_t1_l, &Bltdc_t1_l);
  tree->SetBranchAddress("ltdc_t1_r", &ltdc_t1_r, &Bltdc_t1_r);
  tree->SetBranchAddress("ltdc_bht_l", &ltdc_bht_l, &Bltdc_bht_l);
  tree->SetBranchAddress("ltdc_bht_r", &ltdc_bht_r, &Bltdc_bht_r);

  tree->SetBranchAddress("tot_bdc_l1", &tot_bdc_l1, &Btot_bdc_l1);
  tree->SetBranchAddress("tot_bdc_l2", &tot_bdc_l2, &Btot_bdc_l2);
  tree->SetBranchAddress("tot_bdc_l3", &tot_bdc_l3, &Btot_bdc_l3);
  tree->SetBranchAddress("tot_bdc_l4", &tot_bdc_l4, &Btot_bdc_l4);
  tree->SetBranchAddress("tot_bdc_l5", &tot_bdc_l5, &Btot_bdc_l5);
  tree->SetBranchAddress("tot_bdc_l6", &tot_bdc_l6, &Btot_bdc_l6);
  tree->SetBranchAddress("tot_bdc_l7", &tot_bdc_l7, &Btot_bdc_l7);
  tree->SetBranchAddress("tot_bdc_l8", &tot_bdc_l8, &Btot_bdc_l8);

  tree->SetBranchAddress("ltdc_bdc_l1", &ltdc_bdc_l1, &Bltdc_bdc_l1);
  tree->SetBranchAddress("ltdc_bdc_l2", &ltdc_bdc_l2, &Bltdc_bdc_l2);
  tree->SetBranchAddress("ltdc_bdc_l3", &ltdc_bdc_l3, &Bltdc_bdc_l3);
  tree->SetBranchAddress("ltdc_bdc_l4", &ltdc_bdc_l4, &Bltdc_bdc_l4);
  tree->SetBranchAddress("ltdc_bdc_l5", &ltdc_bdc_l5, &Bltdc_bdc_l5);
  tree->SetBranchAddress("ltdc_bdc_l6", &ltdc_bdc_l6, &Bltdc_bdc_l6);
  tree->SetBranchAddress("ltdc_bdc_l7", &ltdc_bdc_l7, &Bltdc_bdc_l7);
  tree->SetBranchAddress("ltdc_bdc_l8", &ltdc_bdc_l8, &Bltdc_bdc_l8);

  tree->SetBranchAddress("tot_kldc_l1", &tot_kldc_l1, &Btot_kldc_l1);
  tree->SetBranchAddress("tot_kldc_l2", &tot_kldc_l2, &Btot_kldc_l2);
  tree->SetBranchAddress("tot_kldc_l3", &tot_kldc_l3, &Btot_kldc_l3);
  tree->SetBranchAddress("tot_kldc_l4", &tot_kldc_l4, &Btot_kldc_l4);
  tree->SetBranchAddress("tot_kldc_l5", &tot_kldc_l5, &Btot_kldc_l5);
  tree->SetBranchAddress("tot_kldc_l6", &tot_kldc_l6, &Btot_kldc_l6);
  tree->SetBranchAddress("tot_kldc_l7", &tot_kldc_l7, &Btot_kldc_l7);
  tree->SetBranchAddress("tot_kldc_l8", &tot_kldc_l8, &Btot_kldc_l8);

  tree->SetBranchAddress("ltdc_kldc_l1", &ltdc_kldc_l1, &Bltdc_kldc_l1);
  tree->SetBranchAddress("ltdc_kldc_l2", &ltdc_kldc_l2, &Bltdc_kldc_l2);
  tree->SetBranchAddress("ltdc_kldc_l3", &ltdc_kldc_l3, &Bltdc_kldc_l3);
  tree->SetBranchAddress("ltdc_kldc_l4", &ltdc_kldc_l4, &Bltdc_kldc_l4);
  tree->SetBranchAddress("ltdc_kldc_l5", &ltdc_kldc_l5, &Bltdc_kldc_l5);
  tree->SetBranchAddress("ltdc_kldc_l6", &ltdc_kldc_l6, &Bltdc_kldc_l6);
  tree->SetBranchAddress("ltdc_kldc_l7", &ltdc_kldc_l7, &Bltdc_kldc_l7);
  tree->SetBranchAddress("ltdc_kldc_l8", &ltdc_kldc_l8, &Bltdc_kldc_l8);

  tree->SetBranchAddress("tot_bft_l1", &tot_bft_l1, &Btot_bft_l1);
  tree->SetBranchAddress("tot_bft_l2", &tot_bft_l2, &Btot_bft_l2);
  tree->SetBranchAddress("tot_bft_l3", &tot_bft_l3, &Btot_bft_l3);
  tree->SetBranchAddress("tot_bft_l4", &tot_bft_l4, &Btot_bft_l4);
  tree->SetBranchAddress("tot_bft_l5", &tot_bft_l5, &Btot_bft_l5);
  tree->SetBranchAddress("tot_bft_l6", &tot_bft_l6, &Btot_bft_l6);

  tree->SetBranchAddress("ltdc_bft_l1", &ltdc_bft_l1, &Bltdc_bft_l1);
  tree->SetBranchAddress("ltdc_bft_l2", &ltdc_bft_l2, &Bltdc_bft_l2);
  tree->SetBranchAddress("ltdc_bft_l3", &ltdc_bft_l3, &Bltdc_bft_l3);
  tree->SetBranchAddress("ltdc_bft_l4", &ltdc_bft_l4, &Bltdc_bft_l4);
  tree->SetBranchAddress("ltdc_bft_l5", &ltdc_bft_l5, &Bltdc_bft_l5);
  tree->SetBranchAddress("ltdc_bft_l6", &ltdc_bft_l6, &Bltdc_bft_l6);

  tree->SetBranchAddress("tot_sft_l1", &tot_sft_l1, &Btot_sft_l1);
  tree->SetBranchAddress("tot_sft_l2", &tot_sft_l2, &Btot_sft_l2);
  tree->SetBranchAddress("tot_sft_l3", &tot_sft_l3, &Btot_sft_l3);
  tree->SetBranchAddress("tot_sft_l4", &tot_sft_l4, &Btot_sft_l4);
  tree->SetBranchAddress("tot_sft_l5", &tot_sft_l5, &Btot_sft_l5);
  tree->SetBranchAddress("tot_sft_l6", &tot_sft_l6, &Btot_sft_l6);

  tree->SetBranchAddress("ltdc_sft_l1", &ltdc_sft_l1, &Bltdc_sft_l1);
  tree->SetBranchAddress("ltdc_sft_l2", &ltdc_sft_l2, &Bltdc_sft_l2);
  tree->SetBranchAddress("ltdc_sft_l3", &ltdc_sft_l3, &Bltdc_sft_l3);
  tree->SetBranchAddress("ltdc_sft_l4", &ltdc_sft_l4, &Bltdc_sft_l4);
  tree->SetBranchAddress("ltdc_sft_l5", &ltdc_sft_l5, &Bltdc_sft_l5);
  tree->SetBranchAddress("ltdc_sft_l6", &ltdc_sft_l6, &Bltdc_sft_l6);
  
  tree->GetEntry(i_entry);

  for(int i=0; i<NumOfSegUTOF; i++){
     AddHodoRHit(UTOFRHC, DetIdUTOF, 0, i+1,
                 ltdc_utof_l->at(i), tot_utof_l->at(i), 
                 ltdc_utof_r->at(i), tot_utof_r->at(i) );
  }

  for(int i=0; i<NumOfSegDTOF; i++){
     AddHodoRHit(DTOFRHC, DetIdDTOF, 0, i+1,
                 ltdc_dtof_l->at(i), tot_dtof_l->at(i), 
                 ltdc_dtof_r->at(i), tot_dtof_r->at(i) );
  }

  for(int i=0; i<NumOfSegLTOF; i++){
     AddHodoRHit(LTOFRHC, DetIdLTOF, 0, i+1,
                 ltdc_ltof_l->at(i), tot_ltof_l->at(i), 
                 ltdc_ltof_r->at(i), tot_ltof_r->at(i) );
  }

  for(int i=0; i<NumOfSegT0; i++){
     AddHodoRHit(T0RHC, DetIdT0, 0, i+1,
                 ltdc_t0_l->at(i), tot_t0_l->at(i), 
                 ltdc_t0_r->at(i), tot_t0_r->at(i) );
  }

  for(int i=0; i<NumOfSegT0r; i++){
     AddHodoRHit(T0rRHC, DetIdT0r, 0, i+1,
                 ltdc_t0r_l->at(i), tot_t0r_l->at(i), 
                 ltdc_t0r_r->at(i), tot_t0r_r->at(i) );
  }

  for(int i=0; i<NumOfSegBref; i++){
     if(i==0){
        AddHodoRHit(BrefRHC, DetIdBref, 0, i+1,
                    ltdc_bref_l->at(i), tot_bref_l->at(i), 
                    ltdc_bref_l->at(i), tot_bref_l->at(i) );
     }
     if(i==1){
        AddHodoRHit(BrefRHC, DetIdBref, 0, i-1+1,
                    ltdc_bref_r->at(i-1), tot_bref_r->at(i-1), 
                    ltdc_bref_r->at(i-1), tot_bref_r->at(i-1) );
     }
  }
  
  for(int i=0; i<NumOfSegT1; i++){
     AddHodoRHit(T1RHC, DetIdT1, 0, i+1,
                 ltdc_t1_l->at(i), tot_t1_l->at(i), 
                 ltdc_t1_r->at(i), tot_t1_r->at(i) );
  }

  for(int i=0; i<NumOfSegBHT; i++){
     AddHodoRHit(BHTRHC, DetIdBHT, 0, i+1,
                 ltdc_bht_l->at(i), tot_bht_l->at(i), 
                 ltdc_bht_r->at(i), tot_bht_r->at(i) );
  }

  ////BDC
  for(int i=0; i<NumOfWireBDCX; i++)
     AddDCRHit(BDCRHC[1], 1, i+1, ltdc_bdc_l1->at(i), tot_bdc_l1->at(i) );

  for(int i=0; i<NumOfWireBDCX; i++)
     AddDCRHit(BDCRHC[2], 2, i+1, ltdc_bdc_l2->at(i), tot_bdc_l2->at(i) );

  for(int i=0; i<NumOfWireBDCUV; i++)
     AddDCRHit(BDCRHC[3], 3, i+1, ltdc_bdc_l3->at(i), tot_bdc_l3->at(i) );

  for(int i=0; i<NumOfWireBDCUV; i++)
     AddDCRHit(BDCRHC[4], 4, i+1, ltdc_bdc_l4->at(i), tot_bdc_l4->at(i) );

  for(int i=0; i<NumOfWireBDCX; i++)
     AddDCRHit(BDCRHC[5], 5, i+1, ltdc_bdc_l5->at(i), tot_bdc_l5->at(i) );

  for(int i=0; i<NumOfWireBDCX; i++)
     AddDCRHit(BDCRHC[6], 6, i+1, ltdc_bdc_l6->at(i), tot_bdc_l6->at(i) );

  for(int i=0; i<NumOfWireBDCUV; i++)
     AddDCRHit(BDCRHC[7], 7, i+1, ltdc_bdc_l7->at(i), tot_bdc_l7->at(i) );

  for(int i=0; i<NumOfWireBDCUV; i++)
     AddDCRHit(BDCRHC[8], 8, i+1, ltdc_bdc_l8->at(i), tot_bdc_l8->at(i) );

  ////KLDC
  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[1], 1, i+1, ltdc_kldc_l1->at(i), tot_kldc_l1->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[2], 2, i+1, ltdc_kldc_l2->at(i), tot_kldc_l2->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[3], 3, i+1, ltdc_kldc_l3->at(i), tot_kldc_l3->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[4], 4, i+1, ltdc_kldc_l4->at(i), tot_kldc_l4->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[5], 5, i+1, ltdc_kldc_l5->at(i), tot_kldc_l5->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[6], 6, i+1, ltdc_kldc_l6->at(i), tot_kldc_l6->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[7], 7, i+1, ltdc_kldc_l7->at(i), tot_kldc_l7->at(i) );

  for(int i=0; i<NumOfWireKLDC; i++)
     AddDCRHit(KLDCRHC[8], 8, i+1, ltdc_kldc_l8->at(i), tot_kldc_l8->at(i) );

  ////BFT
  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[1], 1, i+1, ltdc_bft_l1->at(i), tot_bft_l1->at(i) );

  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[2], 2, i+1, ltdc_bft_l2->at(i), tot_bft_l2->at(i) );

  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[3], 3, i+1, ltdc_bft_l3->at(i), tot_bft_l3->at(i) );

  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[4], 4, i+1, ltdc_bft_l4->at(i), tot_bft_l4->at(i) );

  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[5], 5, i+1, ltdc_bft_l5->at(i), tot_bft_l5->at(i) );

  // for(int i=0; i<NumOfFiberBFT; i++)
  //    AddTrRHit(BFTRHC[6], 6, i+1, ltdc_bft_l6->at(i), tot_bft_l6->at(i) );

  ////SFT
  // for(int i=0; i<NumOfFiberSFTX; i++)
  //    AddTrRHit(SFTRHC[1], 1, i+1, ltdc_sft_l1->at(i), tot_sft_l1->at(i) );

  // for(int i=0; i<NumOfFiberSFTUV; i++)
  //    AddTrRHit(SFTRHC[2], 2, i+1, ltdc_sft_l2->at(i), tot_sft_l2->at(i) );

  // for(int i=0; i<NumOfFiberSFTUV; i++)
  //    AddTrRHit(SFTRHC[3], 3, i+1, ltdc_sft_l3->at(i), tot_sft_l3->at(i) );

  // for(int i=0; i<NumOfFiberSFTUV; i++)
  //    AddTrRHit(SFTRHC[4], 4, i+1, ltdc_sft_l4->at(i), tot_sft_l4->at(i) );

  // for(int i=0; i<NumOfFiberSFTUV; i++)
  //    AddTrRHit(SFTRHC[5], 5, i+1, ltdc_sft_l5->at(i), tot_sft_l5->at(i) );

  // for(int i=0; i<NumOfFiberSFTX; i++)
  //    AddTrRHit(SFTRHC[6], 6, i+1, ltdc_sft_l6->at(i), tot_sft_l6->at(i) );

  //Clear processe
  delete tot_utof_l;
  delete tot_utof_r;
  delete tot_dtof_l;
  delete tot_dtof_r;
  delete tot_ltof_l;
  delete tot_ltof_r;
  delete tot_t0_l;
  delete tot_t0_r;
  delete tot_t0r_l;
  delete tot_t0r_r;
  delete tot_bref_l;
  delete tot_bref_r;
  delete tot_t1_l;
  delete tot_t1_r;
  delete tot_bht_l;
  delete tot_bht_r;

  delete ltdc_utof_l;
  delete ltdc_utof_r;
  delete ltdc_dtof_l;
  delete ltdc_dtof_r;
  delete ltdc_ltof_l;
  delete ltdc_ltof_r;
  delete ltdc_t0_l;
  delete ltdc_t0_r;
  delete ltdc_t0r_l;
  delete ltdc_t0r_r;
  delete ltdc_bref_l;
  delete ltdc_bref_r;
  delete ltdc_t1_l;
  delete ltdc_t1_r;
  delete ltdc_bht_l;
  delete ltdc_bht_r;

  delete tot_bdc_l1;
  delete tot_bdc_l2;
  delete tot_bdc_l3;
  delete tot_bdc_l4;
  delete tot_bdc_l5;
  delete tot_bdc_l6;
  delete tot_bdc_l7;
  delete tot_bdc_l8;

  delete ltdc_bdc_l1;
  delete ltdc_bdc_l2;
  delete ltdc_bdc_l3;
  delete ltdc_bdc_l4;
  delete ltdc_bdc_l5;
  delete ltdc_bdc_l6;
  delete ltdc_bdc_l7;
  delete ltdc_bdc_l8;

  delete tot_kldc_l1;
  delete tot_kldc_l2;
  delete tot_kldc_l3;
  delete tot_kldc_l4;
  delete tot_kldc_l5;
  delete tot_kldc_l6;
  delete tot_kldc_l7;
  delete tot_kldc_l8;

  delete ltdc_kldc_l1;
  delete ltdc_kldc_l2;
  delete ltdc_kldc_l3;
  delete ltdc_kldc_l4;
  delete ltdc_kldc_l5;
  delete ltdc_kldc_l6;
  delete ltdc_kldc_l7;
  delete ltdc_kldc_l8;

  delete tot_bft_l1;
  delete tot_bft_l2;
  delete tot_bft_l3;
  delete tot_bft_l4;
  delete tot_bft_l5;
  delete tot_bft_l6;

  delete ltdc_bft_l1;
  delete ltdc_bft_l2;
  delete ltdc_bft_l3;
  delete ltdc_bft_l4;
  delete ltdc_bft_l5;
  delete ltdc_bft_l6;

  delete tot_sft_l1;
  delete tot_sft_l2;
  delete tot_sft_l3;
  delete tot_sft_l4;
  delete tot_sft_l5;
  delete tot_sft_l6;

  delete ltdc_sft_l1;
  delete ltdc_sft_l2;
  delete ltdc_sft_l3;
  delete ltdc_sft_l4;
  delete ltdc_sft_l5;
  delete ltdc_sft_l6;
  
  tree->ResetBranchAddress(Btot_utof_l);
  tree->ResetBranchAddress(Btot_utof_r);
  tree->ResetBranchAddress(Btot_dtof_l);
  tree->ResetBranchAddress(Btot_dtof_r);
  tree->ResetBranchAddress(Btot_ltof_l);
  tree->ResetBranchAddress(Btot_ltof_r);
  tree->ResetBranchAddress(Btot_t0_l);
  tree->ResetBranchAddress(Btot_t0_r);
  tree->ResetBranchAddress(Btot_t0r_l);
  tree->ResetBranchAddress(Btot_t0r_r);
  tree->ResetBranchAddress(Btot_bref_l);
  tree->ResetBranchAddress(Btot_bref_r);
  tree->ResetBranchAddress(Btot_t1_l);
  tree->ResetBranchAddress(Btot_t1_r);
  tree->ResetBranchAddress(Btot_bht_l);
  tree->ResetBranchAddress(Btot_bht_r);

  tree->ResetBranchAddress(Bltdc_utof_l);
  tree->ResetBranchAddress(Bltdc_utof_r);
  tree->ResetBranchAddress(Bltdc_dtof_l);
  tree->ResetBranchAddress(Bltdc_dtof_r);
  tree->ResetBranchAddress(Bltdc_ltof_l);
  tree->ResetBranchAddress(Bltdc_ltof_r);
  tree->ResetBranchAddress(Bltdc_t0_l);
  tree->ResetBranchAddress(Bltdc_t0_r);
  tree->ResetBranchAddress(Bltdc_t0r_l);
  tree->ResetBranchAddress(Bltdc_t0r_r);
  tree->ResetBranchAddress(Bltdc_bref_l);
  tree->ResetBranchAddress(Bltdc_bref_r);
  tree->ResetBranchAddress(Bltdc_t1_l);
  tree->ResetBranchAddress(Bltdc_t1_r);
  tree->ResetBranchAddress(Bltdc_bht_l);
  tree->ResetBranchAddress(Bltdc_bht_r);

  tree->ResetBranchAddress(Btot_bdc_l1);
  tree->ResetBranchAddress(Btot_bdc_l2);
  tree->ResetBranchAddress(Btot_bdc_l3);
  tree->ResetBranchAddress(Btot_bdc_l4);
  tree->ResetBranchAddress(Btot_bdc_l5);
  tree->ResetBranchAddress(Btot_bdc_l6);
  tree->ResetBranchAddress(Btot_bdc_l7);
  tree->ResetBranchAddress(Btot_bdc_l8);

  tree->ResetBranchAddress(Bltdc_bdc_l1);
  tree->ResetBranchAddress(Bltdc_bdc_l2);
  tree->ResetBranchAddress(Bltdc_bdc_l3);
  tree->ResetBranchAddress(Bltdc_bdc_l4);
  tree->ResetBranchAddress(Bltdc_bdc_l5);
  tree->ResetBranchAddress(Bltdc_bdc_l6);
  tree->ResetBranchAddress(Bltdc_bdc_l7);
  tree->ResetBranchAddress(Bltdc_bdc_l8);

  tree->ResetBranchAddress(Btot_kldc_l1);
  tree->ResetBranchAddress(Btot_kldc_l2);
  tree->ResetBranchAddress(Btot_kldc_l3);
  tree->ResetBranchAddress(Btot_kldc_l4);
  tree->ResetBranchAddress(Btot_kldc_l5);
  tree->ResetBranchAddress(Btot_kldc_l6);
  tree->ResetBranchAddress(Btot_kldc_l7);
  tree->ResetBranchAddress(Btot_kldc_l8);

  tree->ResetBranchAddress(Bltdc_kldc_l1);
  tree->ResetBranchAddress(Bltdc_kldc_l2);
  tree->ResetBranchAddress(Bltdc_kldc_l3);
  tree->ResetBranchAddress(Bltdc_kldc_l4);
  tree->ResetBranchAddress(Bltdc_kldc_l5);
  tree->ResetBranchAddress(Bltdc_kldc_l6);
  tree->ResetBranchAddress(Bltdc_kldc_l7);
  tree->ResetBranchAddress(Bltdc_kldc_l8);

  tree->ResetBranchAddress(Btot_bft_l1);
  tree->ResetBranchAddress(Btot_bft_l2);
  tree->ResetBranchAddress(Btot_bft_l3);
  tree->ResetBranchAddress(Btot_bft_l4);
  tree->ResetBranchAddress(Btot_bft_l5);
  tree->ResetBranchAddress(Btot_bft_l6);

  tree->ResetBranchAddress(Bltdc_bft_l1);
  tree->ResetBranchAddress(Bltdc_bft_l2);
  tree->ResetBranchAddress(Bltdc_bft_l3);
  tree->ResetBranchAddress(Bltdc_bft_l4);
  tree->ResetBranchAddress(Bltdc_bft_l5);
  tree->ResetBranchAddress(Bltdc_bft_l6);

  tree->ResetBranchAddress(Btot_sft_l1);
  tree->ResetBranchAddress(Btot_sft_l2);
  tree->ResetBranchAddress(Btot_sft_l3);
  tree->ResetBranchAddress(Btot_sft_l4);
  tree->ResetBranchAddress(Btot_sft_l5);
  tree->ResetBranchAddress(Btot_sft_l6);

  tree->ResetBranchAddress(Bltdc_sft_l1);
  tree->ResetBranchAddress(Bltdc_sft_l2);
  tree->ResetBranchAddress(Bltdc_sft_l3);
  tree->ResetBranchAddress(Bltdc_sft_l4);
  tree->ResetBranchAddress(Bltdc_sft_l5);
  tree->ResetBranchAddress(Bltdc_sft_l6);

  return true;
}
// RawData::DecodeRawHits


const DCRHitContainer & RawData::GetBDCRHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBDC ) layer=0;
  return BDCRHC[layer];
}

const DCRHitContainer & RawData::GetKLDCRHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersKLDC ) layer=0;
  return KLDCRHC[layer];
}

const TrRHitContainer & RawData::GetBFTRHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBFT ) layer=0;
  return BFTRHC[layer];
}

const TrRHitContainer & RawData::GetSFTRHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=0;
  return SFTRHC[layer];
}

const HodoRHitContainer& RawData::GetUTOFRHC() const
{
  return UTOFRHC;
}

const HodoRHitContainer& RawData::GetDTOFRHC() const
{
  return DTOFRHC;
}

const HodoRHitContainer& RawData::GetLTOFRHC() const
{
  return LTOFRHC;
}

const HodoRHitContainer& RawData::GetT0RHC() const
{
  return T0RHC;
}

const HodoRHitContainer& RawData::GetT0rRHC() const
{
  return T0rRHC;
}

const HodoRHitContainer& RawData::GetBrefRHC() const
{
  return BrefRHC;
}

const HodoRHitContainer& RawData::GetT1RHC() const
{
  return T1RHC;
}

const HodoRHitContainer& RawData::GetBHTRHC() const
{
  return BHTRHC;
}


