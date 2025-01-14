/*
  DCAnalyzer.cc

  2024/04  K.Shirotori
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

#include "DetectorInfo.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

const double Deg2Rad = acos(-1.)/180. ;
  
DCAnalyzer::DCAnalyzer()
{
}

DCAnalyzer::~DCAnalyzer()
{
}

bool DCAnalyzer::DecodeBDCRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeBDCRawHits]";
  ConfMan *confMan = ConfMan::GetConfManager();
  const double dctdclow = confMan->BDCTRangeLow();
  const double dctdchigh = confMan->BDCTRangeHigh();
  const double dctdctot = confMan->BDCTRangeTOT();

  clearBDCHits();

  for( int layer=1; layer<=NumOfLayersBDC; ++layer ){
     const DCRHitContainer &cont=rawData->GetBDCRHC(layer);
     int nh=cont.size();
     for( int i=0; i<nh; ++i ){
        DCRawHit *rhit=cont[i];
      
        DCHit *hit=new DCHit(rhit->LayerId()+PlOffsBDC, rhit->WireId());

        int nhltdc=rhit->GetSize_lTdc();
        int nhtot=rhit->GetSize_Tot();
        for( int j=0; j<nhltdc; ++j ){
           double ltdc = rhit->GetlTdc(j); 
           double tdc1st = 9999.0;
           if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc < tdc1st ){
              tdc1st=ltdc;
              if( nhtot==nhtot ){
                 //std::cout<< nhltdc << "<-->" << nhtot << std::endl;
                 double tot = rhit->GetTot(j); 
                 if( tot>dctdctot ) hit->SetTdcVal( tdc1st );
              }
           }
        }
        if(!hit) continue;
        
        if(hit->CalcDCObservables())
           BDCHC[layer].push_back(hit);
        else
           delete hit;
     }
  }
  
  return true;
}


bool DCAnalyzer::DecodeKLDCRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeKLDCRawHits]";
  ConfMan *confMan = ConfMan::GetConfManager();
  const double dctdclow = confMan->KLDCTRangeLow();
  const double dctdchigh = confMan->KLDCTRangeHigh();
  const double dctdctot = confMan->KLDCTRangeTOT();

  clearKLDCHits();

  for( int layer=1; layer<=NumOfLayersKLDC; ++layer ){
     const DCRHitContainer &cont=rawData->GetKLDCRHC(layer);
     int nh=cont.size();
     for( int i=0; i<nh; ++i ){
        DCRawHit *rhit=cont[i];
      
        DCHit *hit=new DCHit(rhit->LayerId()+PlOffsKLDC, rhit->WireId());

        int nhltdc=rhit->GetSize_lTdc();
        int nhtot=rhit->GetSize_Tot();
        for( int j=0; j<nhltdc; ++j ){
           double ltdc = rhit->GetlTdc(j); 
           double tdc1st = 9999.0;
           if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc < tdc1st ){
              tdc1st=ltdc;
              if( nhtot==nhtot ){
                 //std::cout<< nhltdc << "<-->" << nhtot << std::endl;
                 double tot = rhit->GetTot(j); 
                 if( tot>dctdctot ) hit->SetTdcVal( tdc1st );
              }
           }
        }
        if(!hit) continue;
        
        if(hit->CalcDCObservables())
           KLDCHC[layer].push_back(hit);
        else
           delete hit;
     }
  }
  
  return true;
}

bool DCAnalyzer::TrackSearchBDC( void )
{
   ConfMan *confMan = ConfMan::GetConfManager();
   const int trmode = confMan->BDCTRMode();

   clearTracksBDC();

   if(trmode==1){
      int ntrack =
         LocalTrackSearch2( &(BDCHC[1]), TrackBDCCol, 
                            NumOfLayersBDC, MinNumOfHitsBDC );
   }
   else{
      int ntrack =
         LocalTrackSearch( BDCHC, PPInfoBDC, NPPInfoBDC,
                           TrackBDCCol, MinNumOfHitsBDC );
   }
   
   return true;
}

bool DCAnalyzer::TrackSearchKLDC( void )
{
   ConfMan *confMan = ConfMan::GetConfManager();
   const int trmode = confMan->KLDCTRMode();
   
   clearTracksKLDC();
   
   if(trmode==1){
      int ntrack =
         LocalTrackSearch2( &(KLDCHC[1]), TrackKLDCCol, 
                            NumOfLayersKLDC, MinNumOfHitsKLDC );
   }
   else{
      int ntrack =
         LocalTrackSearch( KLDCHC, PPInfoKLDC, NPPInfoKLDC,
                           TrackKLDCCol, MinNumOfHitsKLDC );
   }
   
   return true;
}

void DCAnalyzer::clearBDCHits( void )
{
   for( int l=0; l<=NumOfLayersBDC; ++l ){
      std::for_each( BDCHC[l].begin(),  BDCHC[l].end(),  DeleteObject() );
      BDCHC[l].clear();
   }
}

void DCAnalyzer::clearKLDCHits( void )
{
   for( int l=0; l<=NumOfLayersKLDC; ++l ){
      std::for_each( KLDCHC[l].begin(),  KLDCHC[l].end(),  DeleteObject() );
      KLDCHC[l].clear();
   }
}

void DCAnalyzer::clearTracksBDC( void )
{
   std::for_each( TrackBDCCol.begin(), TrackBDCCol.end(), DeleteObject() );
   TrackBDCCol.clear();
}

void DCAnalyzer::clearTracksKLDC( void )
{
   std::for_each( TrackKLDCCol.begin(), TrackKLDCCol.end(), DeleteObject() );
   TrackKLDCCol.clear();
}

bool DCAnalyzer::ReCalcBDCHits( bool applyRecursively )
{
   for( int l=0; l<=NumOfLayersBDC; ++l ){
      int n=BDCHC[l].size();
      for( int i=0; i<n; ++i ){
         DCHit *hit=(BDCHC[l])[i];
         if(hit) hit->ReCalcDC(applyRecursively);
      }
   }
   
   return true;
}

bool DCAnalyzer::ReCalcKLDCHits( bool applyRecursively )
{
   for( int l=0; l<=NumOfLayersKLDC; ++l ){
      int n=KLDCHC[l].size();
      for( int i=0; i<n; ++i ){
         DCHit *hit=(KLDCHC[l])[i];
         if(hit) hit->ReCalcDC(applyRecursively);
      }
   }
   
   return true;
}

bool DCAnalyzer::ReCalcTrackBDC( bool applyRecursively )
{
   int n=TrackBDCCol.size();
   for( int i=0; i<n; ++i ){
      DCLocalTrack *track=TrackBDCCol[i];
      if( track ) track->ReCalc( applyRecursively );
   }
   return true;
}

bool DCAnalyzer::ReCalcTrackKLDC( bool applyRecursively )
{
   int n=TrackKLDCCol.size();
   for( int i=0; i<n; ++i ){
      DCLocalTrack *track=TrackKLDCCol[i];
      if( track ) track->ReCalc( applyRecursively );
   }
   return true;
}

bool DCAnalyzer::ReCalcAll( void )
{
  ReCalcBDCHits();
  ReCalcKLDCHits();

  ReCalcTrackBDC();
  ReCalcTrackKLDC();

  return true;
}

