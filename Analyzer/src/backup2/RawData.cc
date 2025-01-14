/*
  RawData.cc

  2019/08  K.Shirotori
*/

#include "RawData.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "CLHEP/Random/Randomize.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "Decoder.hh"
#include "DetectorID.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
#include "TemplateLib.hh"

#include "ConfMan.hh"
#include "CMapMan.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
  T0RHC(0), DCRHC()
{}

RawData::~RawData()
{
 clearAll();
}

bool RawData::AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int Layer, int Seg, int UorD, 
			   uIntVec lTdc, uIntVec tTdc, uIntVec Width )
{
  static const std::string funcname = "[RawData::AddHodoRHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->LayerId()==Layer &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Layer, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetUorD( UorD );
    p->SetlTdc( lTdc );
    p->SettTdc( tTdc );
    p->SetWidth( Width );

    //std::cout<< "Ch=" << Seg << " :TDC=" << Tdc.size() << std::endl;

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddDCRHit( DCRHitContainer& cont,
			 int Layer, int Wire, 
			 IntVec lTdc, IntVec tTdc )
{
  static const std::string funcname = "[RawData::AddDCRHit]";
  
  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
  	q->WireId()==Wire ){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Layer, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetlTdc( lTdc );
    p->SettTdc( tTdc );

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

  for( int l=0; l<=PlMaxDC; ++l){
    for_each( DCRHC[l].begin(),  DCRHC[l].end(), DeleteObject());
    DCRHC[l].clear();
  }

  return;
}

bool RawData::DecodeRawHits( Decoder& gDec )
{
  clearAll();

  ConfMan *confMan = ConfMan::GetConfManager();
  CMapMan *mapMan=ConfMan::GetConfManager()->GetCMapManager();

  /////////////////////////////////////////////////////////
  //T0 decoding (HR-TDC)
  MZN_HRTDC::Data_t MZNData;
  MZNData = gDec.GetMznHRTdcData();

  std::vector<std::vector<unsigned int> > HRTdcLeading;
  std::vector<std::vector<unsigned int> > HRTdcTrailing;
  std::vector<std::vector<unsigned int> > HRTdcWidth;
  
  HRTdcLeading.resize(NumOfHRTDC);
  HRTdcTrailing.resize(NumOfHRTDC);
  HRTdcWidth.resize(NumOfHRTDC);

  //For DC analysis
  std::vector<std::vector<int> > DCTdcLeading;
  std::vector<std::vector<int> > DCTdcTrailing;
  DCTdcLeading.resize(NumOfHRTDC);
  DCTdcTrailing.resize(NumOfHRTDC);

  //std::cout<< "MZN Event Num = " << MZNData.EventNum << std::endl;

  //std::cout << "*********************" << std::endl;  
  if( MZNData.Nhit>0 || MZNData.tNhit>0 ){
    for( int ich=0; ich<NumOfHRTDC; ich++ ){
      int Caddr = 1, Naddr = 3, Aaddr = -1;
      Aaddr = ich;
      int DetId, PlId, SegId, UorD;
      bool status;

      for( int j=0; j<MZNData.Nhit; j++ ){
	if( ich==MZNData.Ch[j] ){
	  HRTdcLeading[ich].push_back(MZNData.tdc_leading[j]);
	  DCTdcLeading[ich].push_back(MZNData.tdc_leading[j]);
	  //std::cout << MZNData.Ch[j] << " = "  << MZNData.tdc_leading[j] << std::endl;
	}
      }
      for( int j=0; j<MZNData.tNhit; j++ ){
	if( ich==MZNData.tCh[j] ){
	  HRTdcTrailing[ich].push_back(MZNData.tdc_trailing[j]);
	  DCTdcTrailing[ich].push_back(MZNData.tdc_trailing[j]);
	  //std::cout << MZNData.Ch[j] << " = "  << MZNData.tdc_trailing[j] << std::endl;
	}
      }
      
      if( HRTdcLeading[ich].size()>0 ){
	int size_ltdc = HRTdcLeading[ich].size();
	int size_ttdc = HRTdcTrailing[ich].size();
	for( int ilt=0; ilt<size_ltdc; ilt++ ){
	  int ltdc = HRTdcLeading[ich][ilt];
	  for( int itt=0; itt<size_ttdc; itt++ ){
	    int ttdc = HRTdcTrailing[ich][itt];
	    HRTdcWidth[ich].push_back(ltdc-ttdc);
	  }
	}
      }
      
      status = mapMan->GetLogical(Caddr,Naddr,Aaddr,DetId,PlId,SegId,UorD);
      
      if(status){
	AddHodoRHit( T0RHC, DetId, PlId, SegId, UorD, 
		     HRTdcLeading[ich], HRTdcTrailing[ich], HRTdcWidth[ich] );
	
	if( DetId == 3 && PlId<4 ){
	  if(status){
	    AddDCRHit( DCRHC[PlId+1], PlId, SegId, 
		       DCTdcLeading[ich], DCTdcTrailing[ich] );
	  }
	}	    
	
#if 0
	std::cout << "***********************" << std::endl;
	std::cout << "Caddr=" << Caddr
		  << " Naddr=" << Naddr
		  << " Aaddr=" << Aaddr
		  << " --> "
		  << "DetId=" << DetId 
		  << " PlId=" << PlId 
		  << " SegId=" << SegId
		  << " UorD=" << UorD 
		  << std::endl;
#endif
      }
    }
  }
  
  return true;
}

const HodoRHitContainer& RawData::GetT0RHC() const
{
  return T0RHC;
}

const DCRHitContainer & RawData::GetDCRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxDC ) layer=0;
  return DCRHC[layer];
}
