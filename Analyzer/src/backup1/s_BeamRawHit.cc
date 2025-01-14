/*
  s_BeamRawHit.cc

  2016/2
*/

#include "s_BeamRawHit.hh"
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

s_BeamRawHit::~s_BeamRawHit()
{
  clearRegisteredHits();
}

void s_BeamRawHit::clearRegisteredHits( void )
{
  static const std::string funcname="[s_BeamRawHit::clearRegisteredHits]";

  std::for_each(sT0RHC.begin(), sT0RHC.end(), DeleteObject());
  sT0RHC.clear();

  for( int l=0; l<=PlMaxbSSD; ++l){
    for_each( sbSSDRHC[l].begin(),  sbSSDRHC[l].end(), DeleteObject());
    sbSSDRHC[l].clear();
  }
}

bool s_BeamRawHit::SetsTrRHit( int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_BeamRawHit::SetTrRHit]";

  //BFT
  if( Layer>=PlMinbSSD+PlOffsbSSD && Layer<=PlMaxbSSD+PlOffsbSSD ){
    AddsTrRHit(sbSSDRHC[Layer-PlOffsbSSD], Layer, Wire, PosX, PosY, DL);
  }

  return true;
}

bool s_BeamRawHit::SetsHodoRHit( int DetId, int Seg,
				 double TDC0_t, double TDC0_tot,
				 double TDC1_t, double TDC1_tot,
				 double ADC0_t, double ADC0_hgt,
				 double ADC1_t, double ADC1_hgt )
{
  static const std::string funcname="[s_BeamRawHit::SetHodoRHit]";
  
  const TrGeomMan & geomMan=TrGeomMan::GetInstance();
  
  //T0
  // if( DetId == geomMan.GetDetectorId("T0") ){
    AddsHodoRHit(sT0RHC,
		 DetId, Seg, TDC0_t, TDC0_tot, TDC1_t, TDC1_tot,
		 ADC0_t, ADC0_hgt, ADC1_t, ADC1_hgt);
  // }

  return true;
}
bool s_BeamRawHit::AddsTrRHit( s_TrRHitContainer &cont,
			       int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_BeamRawHit::AddTrRHit]";
  
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

bool s_BeamRawHit::AddsHodoRHit( s_HodoRHitContainer &cont,
				 int DetId, int Seg,
				 double TDC0_t, double TDC0_tot,
				 double TDC1_t, double TDC1_tot,
				 double ADC0_t, double ADC0_hgt,
				 double ADC1_t, double ADC1_hgt)
{
  static const std::string funcname="[s_BeamRawHit::AddHodoRHit]";
  
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

const s_HodoRHitContainer& s_BeamRawHit::GetsT0RHC() const
{
  return sT0RHC;
}

const s_TrRHitContainer & s_BeamRawHit::GetsbSSDRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxbSSD ) layer=0;
  return sbSSDRHC[layer];
}
