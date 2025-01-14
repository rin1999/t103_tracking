/*
  HodoAnalyzer.cc

  2024/04 K. Shirotori
*/

#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "HodoCluster.hh"

#include "TemplateLib.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

const double MaxTimeDifT0 = 2.0;

#define Cluster 1

HodoAnalyzer::HodoAnalyzer()
{
}

HodoAnalyzer::~HodoAnalyzer()
{
  clearUTOFHits();
  clearDTOFHits();
  clearLTOFHits();
  clearT0Hits();
  clearT0rHits();
  clearBrefHits();
  clearT1Hits();
  clearBHTHits();
}

void HodoAnalyzer::clearUTOFHits()
{
   std::for_each(UTOFCont.begin(),UTOFCont.end(),DeleteObject());
   UTOFCont.clear();
}

void HodoAnalyzer::clearDTOFHits()
{
   std::for_each(DTOFCont.begin(),DTOFCont.end(),DeleteObject());
   DTOFCont.clear();
}

void HodoAnalyzer::clearLTOFHits()
{
   std::for_each(LTOFCont.begin(),LTOFCont.end(),DeleteObject());
   LTOFCont.clear();
}

void HodoAnalyzer::clearT0Hits()
{
   std::for_each(T0Cont.begin(),T0Cont.end(),DeleteObject());
   T0Cont.clear();
   std::for_each(T0ClCont.begin(),T0ClCont.end(),DeleteObject());
   T0ClCont.clear();
}

void HodoAnalyzer::clearT0rHits()
{
   std::for_each(T0rCont.begin(),T0rCont.end(),DeleteObject());
   T0rCont.clear();
}

void HodoAnalyzer::clearBrefHits()
{
   std::for_each(BrefCont.begin(),BrefCont.end(),DeleteObject());
   BrefCont.clear();
}

void HodoAnalyzer::clearT1Hits( void )
{
   std::for_each( T1Cont.begin(), T1Cont.end(), DeleteObject() );
   T1Cont.clear();
}

void HodoAnalyzer::clearBHTHits( void )
{
   std::for_each( BHTCont.begin(), BHTCont.end(), DeleteObject() );
   BHTCont.clear();
}

bool HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
   DecodeUTOFHits( rawData );
   DecodeDTOFHits( rawData );
   DecodeLTOFHits( rawData );
   DecodeT0Hits( rawData );
   DecodeT0rHits( rawData );
   DecodeBrefHits( rawData );
   DecodeT1Hits( rawData );
   DecodeBHTHits( rawData );
   
   return true;
}

bool HodoAnalyzer::DecodeUTOFHits( RawData *rawData )
{
   clearUTOFHits();

   const HodoRHitContainer &cont=rawData->GetUTOFRHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() ){
          UTOFCont.push_back(hp);
      }
      else{
         delete hp;
      }
   }

   return true;
}

bool HodoAnalyzer::DecodeDTOFHits( RawData *rawData )
{
   clearDTOFHits();
  
   const HodoRHitContainer &cont=rawData->GetDTOFRHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
         DTOFCont.push_back(hp);
      else
         delete hp;
   }
   
   return true;
}

bool HodoAnalyzer::DecodeLTOFHits( RawData *rawData )
{
   clearLTOFHits();
  
   const HodoRHitContainer &cont=rawData->GetLTOFRHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
         LTOFCont.push_back(hp);
      else
         delete hp;
   }
   
   return true;
}

bool HodoAnalyzer::DecodeT0Hits( RawData *rawData )
{
   clearT0Hits();
   
   const HodoRHitContainer &cont=rawData->GetT0RHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
         T0Cont.push_back(hp);
      else
         delete hp;
   }
   
#if Cluster
   MakeUpClusters( T0Cont, T0ClCont, MaxTimeDifT0 );
#endif
   
   return true;
}

bool HodoAnalyzer::DecodeT0rHits( RawData *rawData )
{
   clearT0rHits();
  
   const HodoRHitContainer &cont=rawData->GetT0rRHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
         T0rCont.push_back(hp);
      else
         delete hp;
   }
   
   return true;
}

bool HodoAnalyzer::DecodeBrefHits( RawData *rawData )
{
  clearBrefHits();
  
  const HodoRHitContainer &cont=rawData->GetBrefRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      BrefCont.push_back(hp);
    else
      delete hp;
  }
    
  return true;
}

bool HodoAnalyzer::DecodeT1Hits( RawData *rawData )
{
   clearT1Hits();
  
   const HodoRHitContainer &cont=rawData->GetT1RHC();
   int nh=cont.size();
   for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
         T1Cont.push_back(hp);
      else
         delete hp;
   }
   
   return true;
}

bool HodoAnalyzer::DecodeBHTHits( RawData *rawData )
{
  clearBHTHits();
  
  const HodoRHitContainer &cont=rawData->GetBHTRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      BHTCont.push_back(hp);
    else
      delete hp;
  }
    
  return true;
}

int HodoAnalyzer::
MakeUpClusters( const Hodo2HitContainer & HitCont,
		HodoClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";
  
  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
  
  int nh=HitCont.size();
  
  std::vector <int> flag(nh,0);
  
  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    Hodo2Hit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->CMeanTime();
    Hodo2Hit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      Hodo2Hit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->CMeanTime();
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
        hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      Hodo2Hit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
        if( j==iB ) continue;
        Hodo2Hit *hit=HitCont[j];
        int seg=hit->SegmentId();
        double cmt=hit->CMeanTime();
        if( (abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
            (abs(seg-segB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
          hitC=hit; ++flag[j]; break;
        }
      }
      if(hitC){
        HodoCluster *cluster=new HodoCluster(hitA,hitB,hitC);
        if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
        HodoCluster *cluster=new HodoCluster(hitA,hitB);
        if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      HodoCluster *cluster=new HodoCluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}				  

bool HodoAnalyzer::ReCalcUTOFHits( bool applyRecursively )
{
  int n=UTOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=UTOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcDTOFHits( bool applyRecursively )
{
  int n=DTOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=DTOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcLTOFHits( bool applyRecursively )
{
  int n=LTOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=LTOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcT0Hits( bool applyRecursively )
{
  int n=T0Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=T0Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcT0rHits( bool applyRecursively )
{
  int n=T0rCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=T0rCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBrefHits( bool applyRecursively )
{
  int n=BrefCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=BrefCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcT1Hits( bool applyRecursively )
{
  int n=T1Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=T1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBHTHits( bool applyRecursively )
{
  int n=BHTCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=BHTCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcT0Clusters( bool applyRecursively )
{
  int n=T0ClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=T0ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcAll( void )
{
  ReCalcUTOFHits();
  ReCalcDTOFHits();
  ReCalcLTOFHits();
  ReCalcT0Hits();
  ReCalcT0rHits();
  ReCalcBrefHits();
  ReCalcT1Hits();
  ReCalcBHTHits();
  
  ReCalcT0Clusters();

  return true;
}

