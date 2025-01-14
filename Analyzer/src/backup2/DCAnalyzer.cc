/*
  DCAnalyzer.cc

  2018/12  K.Shirotori
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "TemplateLib.hh"
#include "DCTrackSearch.hh"
#include "ConfMan.hh"

#include "DetectorID.hh"
#include "DCParameters.hh"

const double Deg2Rad = acos(-1.)/180. ;
  
DCAnalyzer::DCAnalyzer()
{
}

DCAnalyzer::~DCAnalyzer()
{
  clearTracksTestDC();
  clearTestDCHits();
}

bool DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeRawHits]";
  ConfMan *confMan = ConfMan::GetConfManager();
  const int dctdclow = confMan->DCTRangeLow();
  const int dctdchigh = confMan->DCTRangeHigh();
  const int dcwidthcutl = confMan->DCWidthCutLow();
  const int dcwidthcuth = confMan->DCWidthCutHigh();

  clearTestDCHits();
  int id=0;
  for( int layer=1; layer<=NumOfLayersDC; ++layer ){
    const DCRHitContainer &cont=rawData->GetDCRHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];
      
      DCHit *hit=new DCHit(rhit->LayerId()+1+PlOffsDC, rhit->WireId()+1);

      bool widthflag=false;
      int nhltdc=rhit->GetSize_lTdc();
      int nhttdc=rhit->GetSize_tTdc();

      int layerId = rhit->LayerId()+1;
      int wireId = rhit->WireId()+1;

      id=MaxWire*(layerId-1)+(wireId-1);
      int id2=id+1;
      // if(layer==1) id2=rhit->WireId();
      // if(layer==2) id2=rhit->WireId()-1;
      // if(layer==3) id2=rhit->WireId()-1;

      for( int irt=0; irt<nhltdc; ++irt ){
	int ltdc=rhit->GetlTdc(irt)*Tdc2Time2;
	for( int itt=0; itt<nhttdc; ++itt ){
	  int ttdc= rhit->GettTdc(itt)*Tdc2Time2;

	  int twidth = -1;
	  if(id2%2==1) twidth = ltdc-ttdc;
	  if(id2%2==0) twidth = ttdc-ltdc;
	  //if(id2%2==1 && id !=15) twidth = ltdc-ttdc;
	  //if(id2%2==0 && id !=16) twidth = ttdc-ltdc;
	  //if(id==16) twidth = ltdc-ttdc;
	  //if(id==15) twidth = ttdc-ltdc;
	  //std::cout<< twidth << std::endl;
	  if( twidth>dcwidthcutl && twidth<dcwidthcuth) widthflag=true;
	}
      }
      //std::cout<< "->" << widthflag << std::endl;

      if(id2%2==1){// && id !=15) {
	for( int j=0; j<nhltdc; ++j ){
	  int ltdc = rhit->GetlTdc(j)*Tdc2Time2; 
	  int tdc1st = -1;
	  if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	  if( widthflag ) hit->SetTdcVal( tdc1st );
	  //std::cout<< rhit->GetlTdc(j)*Tdc2Time2 << std::endl;
	}
	if(!hit) continue;
      }
      if(id2%2==0){// && id != 16) {
	for( int j=0; j<nhttdc; ++j ){
	  int ttdc = rhit->GettTdc(j)*Tdc2Time2; 
	  int tdc1st = -1;
	  if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	  if( widthflag ) hit->SetTdcVal( tdc1st );
	  //std::cout<< rhit->GettTdc(j)*Tdc2Time2 << std::endl;
	}
	if(!hit) continue;
      }
      /*
      if(id==16) {
	for( int j=0; j<nhltdc; ++j ){
	  int ltdc = rhit->GetlTdc(j)*Tdc2Time2; 
	  int tdc1st = -1;
	  if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	  if( widthflag ) hit->SetTdcVal( tdc1st );
	  //std::cout<< rhit->GetlTdc(j)*Tdc2Time2 << std::endl;
	}
	if(!hit) continue;
      }
      if(id== 15) {
	for( int j=0; j<nhttdc; ++j ){
	  int ttdc = rhit->GettTdc(j)*Tdc2Time2; 
	  int tdc1st = -1;
	  if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	  if( widthflag ) hit->SetTdcVal( tdc1st );
	  //std::cout<< rhit->GettTdc(j)*Tdc2Time2 << std::endl;
	}
	if(!hit) continue;
      }
      */




      if(hit->CalcDCObservables())
  	TestDCHC[layer].push_back(hit);
      else
  	delete hit;
    }
  }
      
  return true;
}


bool DCAnalyzer::TrackSearchTestDC( void )
{  
  int ntrack =
    LocalTrackSearch( &(TestDCHC[1]), TrackTestDCCol, 
		      NumOfLayersDC, MinNumOfHitsDC );

  return true;
}

void DCAnalyzer::clearTestDCHits( void )
{
  for( int l=0; l<=NumOfLayersDC; ++l ){
    for_each( TestDCHC[l].begin(),  TestDCHC[l].end(),  DeleteObject() );
    TestDCHC[l].clear();
  }
}

void DCAnalyzer::clearTracksTestDC( void )
{
  for_each( TrackTestDCCol.begin(), TrackTestDCCol.end(), DeleteObject() );
  TrackTestDCCol.clear();
}

bool DCAnalyzer::ReCalcTestDCHits( bool applyRecursively )
{
  for( int l=0; l<=NumOfLayersDC; ++l ){
    int n=TestDCHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(TestDCHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }

  return true;
}

bool DCAnalyzer::ReCalcTrackTestDC( bool applyRecursively )
{
  int n=TrackTestDCCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackTestDCCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}



bool DCAnalyzer::ReCalcAll( void )
{
  ReCalcTestDCHits();

  ReCalcTrackTestDC();

  return true;
}

