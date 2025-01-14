/*
  s_ScatRawHit.cc

  2016/2
*/

#include "s_ScatRawHit.hh"
#include "s_TrRawHit.hh"
#include "s_HodoRawHit.hh"
#include "ConfMan.hh"
#include "TrGeomMan.hh"
#include "DetectorID.hh"
#include "TemplateLib.hh"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>

s_ScatRawHit::~s_ScatRawHit()
{
  clearRegisteredHits();
}

void s_ScatRawHit::clearRegisteredHits( void )
{
  static const std::string funcname="[s_ScatRawHit::clearRegisteredHits]";

  std::for_each(sRPCRHC.begin(), sRPCRHC.end(), DeleteObject());
  sRPCRHC.clear();

  for( int l=0; l<=PlMaxsSSD1; ++l){
    for_each( ssSSD1RHC[l].begin(),  ssSSD1RHC[l].end(), DeleteObject());
    ssSSD1RHC[l].clear();
  }

  for( int l=0; l<=PlMaxsSSD2; ++l){
    for_each( ssSSD2RHC[l].begin(),  ssSSD2RHC[l].end(), DeleteObject());
    ssSSD2RHC[l].clear();
  }

}

bool s_ScatRawHit::SetsTrRHit( int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_ScatRawHit::SetTrRHit]";
  
  //sSSD1
  if( Layer>=PlMinsSSD1+PlOffssSSD1 && Layer<=PlMaxsSSD1+PlOffssSSD1 ){
    AddsTrRHit(ssSSD1RHC[Layer-PlOffssSSD1], Layer, Wire, PosX, PosY, DL);
  }
  
  //sSSD2
  if( Layer>=PlMinsSSD2+PlOffssSSD2 && Layer<=PlMaxsSSD2+PlOffssSSD2 ){
    AddsTrRHit(ssSSD2RHC[Layer-PlOffssSSD2], Layer, Wire, PosX, PosY, DL);
  }
  
  return true;
}

bool s_ScatRawHit::SetsHodoRHit( int DetId, int Seg,
				 double TDC0_t, double TDC0_tot,
				 double TDC1_t, double TDC1_tot,
				 double ADC0_t, double ADC0_hgt,
				 double ADC1_t, double ADC1_hgt )
{
  static const std::string funcname="[s_ScatRawHit::SetHodoRHit]";
  
  const TrGeomMan & geomMan=TrGeomMan::GetInstance();
  
  //RPC
  // if( DetId == geomMan.GetDetectorId("RPC") ){
    AddsHodoRHit(sRPCRHC,
		 DetId, Seg, TDC0_t, TDC0_tot, TDC1_t, TDC1_tot,
		 ADC0_t, ADC0_hgt, ADC1_t, ADC1_hgt);
  // }

  return true;
}

bool s_ScatRawHit::AddsTrRHit( s_TrRHitContainer &cont,
			       int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_ScatRawHit::AddTrRHit]";
  
  s_TrRawHit *p=0;
  int nh=cont.size();
  
  for( int i=0; i<nh; ++i ){
    s_TrRawHit *q=cont[i];
    if( q->LayerId()==Layer && 
	q->WireId()==Wire ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_TrRawHit( Layer, Wire );
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


bool s_ScatRawHit::AddsHodoRHit( s_HodoRHitContainer &cont,
				 int DetId, int Seg,
				 double TDC0_t, double TDC0_tot,
				 double TDC1_t, double TDC1_tot,
				 double ADC0_t, double ADC0_hgt,
				 double ADC1_t, double ADC1_hgt)
{
  static const std::string funcname="[s_ScatRawHit::AddHodoRHit]";
  
  s_HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_HodoRawHit( DetId, Seg );
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

const s_HodoRHitContainer& s_ScatRawHit::GetsRPCRHC() const
{
  return sRPCRHC;
}

const s_TrRHitContainer & s_ScatRawHit::GetssSSD1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxsSSD1 ) layer=0;
  return ssSSD1RHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetssSSD2RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxsSSD2 ) layer=0;
  return ssSSD2RHC[layer];
}
