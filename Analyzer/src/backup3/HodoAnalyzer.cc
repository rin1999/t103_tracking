/*
  HodoAnalyzer.cc
*/

#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "BH2Hit.hh"
#include "HodoCluster.hh"
#include "BH2Cluster.hh"

#include "TemplateLib.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

const double MaxTimeDifBh1 = 2.0;
const double MaxTimeDifBh2 = 2.0;
const double MaxTimeDifTof = 3.5;
const double MaxTimeDifLc  = 5.5;

#define Cluster 1

HodoAnalyzer::HodoAnalyzer()
{
}

HodoAnalyzer::~HodoAnalyzer()
{
  clearACHits();
  clearLCHits();
  clearTOFHits();
  clearBACHits();
  clearBH2Hits();
  clearBH1Hits();
  clearGCHits();
}

void HodoAnalyzer::clearGCHits()
{
  for_each(GCCont.begin(),GCCont.end(),DeleteObject());
}

void HodoAnalyzer::clearBH1Hits()
{
  for_each(BH1Cont.begin(),BH1Cont.end(),DeleteObject());
  for_each(BH1ClCont.begin(),BH1ClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearBH2Hits()
{
  for_each(BH2Cont.begin(),BH2Cont.end(),DeleteObject());
  for_each(BH2ClCont.begin(),BH2ClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearBACHits()
{
  for_each(BACCont.begin(),BACCont.end(),DeleteObject());
}

void HodoAnalyzer::clearTOFHits( void )
{
  for_each( TOFCont.begin(), TOFCont.end(), DeleteObject() );
  for_each( TOFClCont.begin(), TOFClCont.end(), DeleteObject() );
}

void HodoAnalyzer::clearLCHits( void )
{
  for_each( LCCont.begin(), LCCont.end(), DeleteObject() );
  for_each( LCClCont.begin(), LCClCont.end(), DeleteObject() );
}

void HodoAnalyzer::clearACHits( void )
{
  for( int l=0; l<=NumOfLayersAc; ++l )
    for_each( ACCont[l].begin(), ACCont[l].end(), DeleteObject() );
}

bool HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeGCHits( rawData );
  DecodeBH1Hits( rawData );
  DecodeBH2Hits( rawData );
  DecodeBACHits( rawData );
  DecodeTOFHits( rawData );
  DecodeLCHits( rawData );
  DecodeACHits( rawData );

  return true;
}

bool HodoAnalyzer::DecodeGCHits( RawData *rawData )
{
  clearGCHits();
  
  const HodoRHitContainer &cont=rawData->GetGCRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      GCCont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

bool HodoAnalyzer::DecodeBH1Hits( RawData *rawData )
{
  clearBH1Hits();
  
  const HodoRHitContainer &cont=rawData->GetBH1RawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	BH1Cont.push_back(hp);
      else
	delete hp;
    }
  }
  
#if Cluster
  MakeUpClusters( BH1Cont, BH1ClCont, MaxTimeDifBh1 );
#endif
  
  return true;
}

bool HodoAnalyzer::DecodeBH2Hits( RawData *rawData )
{
  clearBH2Hits();
  
  const HodoRHitContainer &cont=rawData->GetBH2RawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      BH2Hit *hp=new BH2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	BH2Cont.push_back(hp);
      else
	delete hp;
    }
  }

  
#if Cluster
  MakeUpClusters( BH2Cont, BH2ClCont, MaxTimeDifBh2 );
#endif
    
  return true;
}

bool HodoAnalyzer::DecodeBACHits( RawData *rawData )
{
  clearBACHits();
  
  const HodoRHitContainer &cont=rawData->GetBACRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    Hodo1Hit *hp=new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      BACCont.push_back(hp);
    else
      delete hp;
  }
    
  return true;
}

bool HodoAnalyzer::DecodeTOFHits( RawData *rawData )
{
  clearTOFHits();

  const HodoRHitContainer &cont=rawData->GetTOFRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	TOFCont.push_back(hp);
      else
	delete hp;
    }
  }

#if Cluster
  MakeUpClusters( TOFCont, TOFClCont, MaxTimeDifTof );
#endif

  return true;
}

bool HodoAnalyzer::DecodeLCHits( RawData *rawData )
{
  clearLCHits();

  const HodoRHitContainer &cont=rawData->GetLCRawHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    int Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
    if( Tu>0 && Td>0 ){
      Hodo2Hit *hp=new Hodo2Hit( hit );
      if( !hp ) continue;
      if( hp->calculate() )
	LCCont.push_back(hp);
      else
	delete hp;
    }
  }

#if Cluster
  MakeUpClusters( LCCont, LCClCont, MaxTimeDifLc );
#endif

  return true;
}

bool HodoAnalyzer::DecodeACHits( RawData *rawData )
{
  static const std::string funcname =
    "HodoAnalyzer::DecodeAcHits";

  clearACHits();

  for( int l=1; l<=NumOfLayersAc; ++l ){ 
    const HodoRHitContainer &cont=rawData->GetACRawHC(l);
    int nh=cont.size();
     for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      if( !hit ) continue;
      int T=hit->GetTdc1();
      if( T>0 ){
	Hodo1Hit *hp=new Hodo1Hit( hit );
	if( !hp ) continue;
	if( hp->calculate() )
	  ACCont[l].push_back(hp);
	else
	  delete hp;
      }
    }
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

int HodoAnalyzer::
MakeUpClusters( const BH2HitContainer & HitCont,
		BH2ClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";

  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );

  int nh=HitCont.size();

  std::vector <int> flag(nh,0);

  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    BH2Hit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->CMeanTime();
    BH2Hit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      BH2Hit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->CMeanTime();
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
	hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      BH2Hit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
        if( j==iB ) continue;
        BH2Hit *hit=HitCont[j];
        int seg=hit->SegmentId();
        double cmt=hit->CMeanTime();
        if( (abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
            (abs(seg-segB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
          hitC=hit; ++flag[j]; break;
        }
      }
      if(hitC){
        BH2Cluster *cluster=new BH2Cluster(hitA,hitB,hitC);
        if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
	BH2Cluster *cluster=new BH2Cluster(hitA,hitB);
	if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      BH2Cluster *cluster=new BH2Cluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}				  

bool HodoAnalyzer::ReCalcGCHits( bool applyRecursively )
{
  int n=GCCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=GCCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH1Hits( bool applyRecursively )
{
  int n=BH1Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=BH1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH2Hits( bool applyRecursively )
{
  int n=BH2Cont.size();
  for( int i=0; i<n; ++i ){
    BH2Hit *hit=BH2Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBACHits( bool applyRecursively )
{
  int n=BACCont.size();
  for( int i=0; i<n; ++i ){
    Hodo1Hit *hit=BACCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcTOFHits( bool applyRecursively )
{
  int n=TOFCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=TOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcLCHits( bool applyRecursively )
{
  int n=LCCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit=LCCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcACHits( bool applyRecursively )
{
  for( int l=0; l<=NumOfLayersAc; ++l ){
    int n=(ACCont[l]).size();
    for( int i=0; i<n; ++i ){
      Hodo1Hit *hit=(ACCont[l])[i];
      if(hit) hit->ReCalc(applyRecursively);
    }
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH1Clusters( bool applyRecursively )
{
  int n=BH1ClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=BH1ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcBH2Clusters( bool applyRecursively )
{
  int n=BH2ClCont.size();
  for( int i=0; i<n; ++i ){
    BH2Cluster *cl=BH2ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcTOFClusters( bool applyRecursively )
{
  int n=TOFClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=TOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcLCClusters( bool applyRecursively )
{
  int n=LCClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=LCClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcAll( void )
{
  ReCalcGCHits();
  ReCalcBH1Hits();
  ReCalcBH2Hits();
  ReCalcBACHits();
  ReCalcTOFHits();
  ReCalcLCHits();
  ReCalcACHits();
  
  ReCalcBH1Clusters();
  ReCalcBH2Clusters();
  ReCalcTOFClusters();
  ReCalcLCClusters();

  return true;
}

