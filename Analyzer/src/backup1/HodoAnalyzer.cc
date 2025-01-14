/*
  HodoAnalyzer.cc

  2016/2  K.Shirotori
*/

#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "HodoHit.hh"
#include "HodoCluster.hh"

#include "TemplateLib.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

const double MaxTimeDifT0  = 0.1;
const double MaxTimeDifTOF = 0.1;
const double MaxTimeDifITOF= 0.1;
const double MaxTimeDifPAD = 0.1;

#define Cluster 1

HodoAnalyzer::HodoAnalyzer()
{
}

HodoAnalyzer::~HodoAnalyzer()
{
  clearVDHits();
  clearMFHits();
  clearPID2Hits();
  clearPID1Hits();
  clearRICHHits();
  clearPADHits();
  clearITOFHits();
  clearTOFHits();
  clearT0Hits();
}

void HodoAnalyzer::clearT0Hits()
{
  for_each(T0Cont.begin(),T0Cont.end(),DeleteObject());
  for_each(T0ClCont.begin(),T0ClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearTOFHits()
{
  for_each(TOFCont.begin(),TOFCont.end(),DeleteObject());
  for_each(TOFClCont.begin(),TOFClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearITOFHits()
{
  for_each(ITOFCont.begin(),ITOFCont.end(),DeleteObject());
  for_each(ITOFClCont.begin(),ITOFClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearPADHits()
{
  for_each(PADCont.begin(),PADCont.end(),DeleteObject());
  for_each(PADClCont.begin(),PADClCont.end(),DeleteObject());
}

void HodoAnalyzer::clearRICHHits()
{
  for_each(RICHCont.begin(),RICHCont.end(),DeleteObject());
}

void HodoAnalyzer::clearPID1Hits()
{
  for_each(PID1Cont.begin(),PID1Cont.end(),DeleteObject());
}

void HodoAnalyzer::clearPID2Hits()
{
  for_each(PID2Cont.begin(),PID2Cont.end(),DeleteObject());
}

void HodoAnalyzer::clearMFHits()
{
  for_each(MFCont.begin(),MFCont.end(),DeleteObject());
}

void HodoAnalyzer::clearVDHits()
{
  for_each(VDCont.begin(),VDCont.end(),DeleteObject());
}

bool HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeT0Hits( rawData );
  DecodeTOFHits( rawData );
  DecodeITOFHits( rawData );
  DecodePADHits( rawData );
  DecodeRICHHits( rawData );
  DecodePID1Hits( rawData );
  DecodePID2Hits( rawData );
  DecodeMFHits( rawData );
  DecodeVDHits( rawData );

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
    HodoHit *hp=new HodoHit( hit );
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

bool HodoAnalyzer::DecodeTOFHits( RawData *rawData )
{
  clearTOFHits();
  
  const HodoRHitContainer &cont=rawData->GetTOFRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      TOFCont.push_back(hp);
    else
      delete hp;
  }
  
#if Cluster
  MakeUpClusters( TOFCont, TOFClCont, MaxTimeDifTOF );
#endif
  
  return true;
}

bool HodoAnalyzer::DecodeITOFHits( RawData *rawData )
{
  clearITOFHits();
  
  const HodoRHitContainer &cont=rawData->GetITOFRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      ITOFCont.push_back(hp);
    else
      delete hp;
  }
  
#if Cluster
  MakeUpClusters( ITOFCont, ITOFClCont, MaxTimeDifITOF );
#endif
  
  return true;
}

bool HodoAnalyzer::DecodePADHits( RawData *rawData )
{
  clearPADHits();
  
  const HodoRHitContainer &cont=rawData->GetPADRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      PADCont.push_back(hp);
    else
      delete hp;
  }
  
#if Cluster
  MakeUpClusters( PADCont, PADClCont, MaxTimeDifPAD );
#endif
  
  return true;
}

bool HodoAnalyzer::DecodeRICHHits( RawData *rawData )
{
  clearRICHHits();
  
  const HodoRHitContainer &cont=rawData->GetRICHRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      RICHCont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

bool HodoAnalyzer::DecodePID1Hits( RawData *rawData )
{
  clearPID1Hits();
  
  const HodoRHitContainer &cont=rawData->GetPID1RHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      PID1Cont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

bool HodoAnalyzer::DecodePID2Hits( RawData *rawData )
{
  clearPID2Hits();
  
  const HodoRHitContainer &cont=rawData->GetPID2RHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      PID2Cont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

bool HodoAnalyzer::DecodeMFHits( RawData *rawData )
{
  clearMFHits();
  
  const HodoRHitContainer &cont=rawData->GetMFRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      MFCont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

bool HodoAnalyzer::DecodeVDHits( RawData *rawData )
{
  clearVDHits();
  
  const HodoRHitContainer &cont=rawData->GetVDRHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      VDCont.push_back(hp);
    else
      delete hp;
  }
  
  return true;
}

int HodoAnalyzer::
MakeUpClusters( const HodoHitContainer & HitCont,
		HodoClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";
  
  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
  
  int nh=HitCont.size();
  
  std::vector <int> flag(nh,0);
  
  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    HodoHit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->GetTime(0);
    HodoHit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      HodoHit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->GetTime(0);
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
	hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      HodoHit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
	if( j==iB ) continue;
	HodoHit *hit=HitCont[j];
	int seg=hit->SegmentId();
	double cmt=hit->GetTime(0);
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

bool HodoAnalyzer::ReCalcT0Hits( bool applyRecursively )
{
  int n=T0Cont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=T0Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcTOFHits( bool applyRecursively )
{
  int n=TOFCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=TOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcITOFHits( bool applyRecursively )
{
  int n=ITOFCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=ITOFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcPADHits( bool applyRecursively )
{
  int n=PADCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=PADCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcRICHHits( bool applyRecursively )
{
  int n=RICHCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=RICHCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcPID1Hits( bool applyRecursively )
{
  int n=PID1Cont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=PID1Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcPID2Hits( bool applyRecursively )
{
  int n=PID2Cont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=PID2Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcMFHits( bool applyRecursively )
{
  int n=MFCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=MFCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcVDHits( bool applyRecursively )
{
  int n=VDCont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=VDCont[i];
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

bool HodoAnalyzer::ReCalcTOFClusters( bool applyRecursively )
{
  int n=TOFClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=TOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcITOFClusters( bool applyRecursively )
{
  int n=ITOFClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=ITOFClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcPADClusters( bool applyRecursively )
{
  int n=PADClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=PADClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}


bool HodoAnalyzer::ReCalcAll( void )
{
  ReCalcT0Hits();
  ReCalcTOFHits();
  ReCalcITOFHits();
  ReCalcPADHits();
  ReCalcRICHHits();
  ReCalcPID1Hits();
  ReCalcPID1Hits();
  ReCalcMFHits();
  ReCalcVDHits();

  ReCalcT0Clusters();  
  ReCalcTOFClusters();
  ReCalcITOFClusters();
  ReCalcPADClusters();

  return true;
}

