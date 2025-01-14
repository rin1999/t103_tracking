/*
  TrAnalyzer.cc

  2024/11  K.Shirotori
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

#include "RawData.hh"
#include "TrRawHit.hh"

#include "TrAnalyzer.hh"
#include "TrHit.hh"
#include "TrCluster.hh"
#include "TrLocalTrack.hh"
#include "TemplateLib.hh"
#include "TrTrackSearch.hh"
#include "ConfMan.hh"

#include "TrParameters.hh"

const double Deg2Rad = acos(-1.)/180. ;

const double MaxTimeDifTr = 10.0;

TrAnalyzer::TrAnalyzer()
{
}

TrAnalyzer::~TrAnalyzer()
{
  clearTracksBFT();
  clearBFTCHits();
  clearBFTHits();
}

bool TrAnalyzer::DecodeBFTRawHits( RawData *rawData )
{
  const std::string funcname = "[TrAnalyzer::DecodeRawHits]";

  ConfMan *confMan = ConfMan::GetConfManager();
  const int tdclow = confMan->BFTTRangeLow();
  const int tdchigh = confMan->BFTTRangeHigh();

  clearBFTHits();
  clearBFTCHits();

  // std::cout<< "*************************** 1" << std::endl;
  for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
    const TrRHitContainer &cont=rawData->GetBFTRHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      TrRawHit *rhit=cont[i];

      int nhtdc= rhit->GetSize_lTdc();
      for( int j=0; j<nhtdc; ++j ){
	int ltdc = rhit->GetlTdc(j);
  int tot  = rhit->GetTot(j);

	if( tdclow<ltdc && ltdc<tdchigh ){
	  TrHit *hit=new TrHit(rhit->LayerId()+PlOffsBFT, rhit->FiberId());      
	  hit->SetTdcVal( ltdc );
    hit->SetTotVal( tot );
	  
	  // std::cout<< "**** 2" << std::endl;
	  // std::cout<< rhit->LayerId()+PlOffsBFT << " "
	  //        << rhit->FiberId() << " "
	  //        << rhit->GetlTdc(j) << std::endl;

	  if(hit->CalcTrObservables()){
	    BFTHC[layer].push_back(hit);
	  }
	  else
	    delete hit;
	}
      }
    }

    // std::cout<< "***** 3 Cluster" << std::endl;
    int cnh = MakeUpClusters( BFTHC[layer], BFTClCont, MaxTimeDifTr );
    // std::cout << "ClusterSize = " << cnh  << std::endl;
    // std::cout<< "**** 4" << std::endl;

    for( int j=0; j<cnh; ++j ){
      TrCluster *chit = GetClusterBFT(j);
      int cs = chit->ClusterSize();
      double mf = chit->CMeanFiber();
      double mt = chit->CMeanTime();
      // std::cout << cs <<" M_Fiber:M_Time = " << mf  << ":" << mt << std::endl;
      
      TrHit *hit=new TrHit();
      hit->SetLayer(layer+PlOffsBFT);
      hit->SetClusterSize(cs);
      hit->SetMeanFiber(mf);
      hit->SetMeanTime(mt);
      hit->SetTdcVal(0.0);
	
	if(hit->ReCalcTrObservables()){
	  BFTCHC[layer].push_back(hit);
	}
	else
	  delete hit;
    }
    
    // std::cout << "****************** 5" << std::endl;
    int nhc = BFTCHC[layer].size();
    for( int j=0; j<nhc; ++j ){
      TrHit *hit = BFTCHC[layer][j];
      // std::cout << "Layer = " << hit->GetLayer() 
      //   	  << " MeanFiber = " << hit->GetMeanFiber() 
      //   	<< " MeanTime = " << hit->GetMeanTime(0) 
      //   	<< std::endl;
    }
  }
 
  return true;
}

int TrAnalyzer::MakeUpClusters( const TrHitContainer & HitCont, TrClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "TrAnalyzer::MakeUpClusters";

  if( !ClusterCont.empty() ){
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
    ClusterCont.clear();
  }

  int nh=HitCont.size();
  
  std::vector <int> flag(nh,0);

  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    TrHit *hitA=HitCont[i];
    int fiberA=hitA->GetFiber();
    double cmtA=hitA->GetTime(0); //Only one hits case

    TrHit *hitB=0;
    int iB=-1;
    int fiberB;
    double cmtB;

    for( int j=i+1; j<nh; ++j ){
      TrHit *hit=HitCont[j];
      int fiber=hit->GetFiber(); 
      double cmt=hit->GetTime(0); //Only one hits case
      // std::cout<< "**2" << std::endl;
      // std::cout<< "LayerA = " << hitA->GetLayer() << std::endl;
      // std::cout<< fiber << " - " << fiberA << " = " << abs(fiber-fiberA) << std::endl;
      if( abs(fiber-fiberA)==1 ){
	//std::cout<< "***3" << std::endl;
	hitB=hit; ++flag[j]; iB=j; fiberB=fiber; break;
      }
      // if( abs(fiber-fiberA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
      // 	//std::cout<< "***3" << std::endl;
      // 	hitB=hit; ++flag[j]; iB=j; fiberB=fiber; cmtB=cmt; break;
      // }
    }
    if(hitB){
      TrHit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
	if( j==iB ) continue;
	TrHit *hit=HitCont[j];
	int fiber=hit->GetFiber();
	double cmt=hit->GetTime(0); //Only one hits case
	//std::cout<< "***3 LayerA = " << hitA->GetLayer() << std::endl;
	//std::cout<< "***3 LayerB = " << hitB->GetLayer() << std::endl;
	//std::cout<< fiber << " - " << fiberA << " = " << abs(fiber-fiberA) << std::endl;
	//std::cout<< fiber << " - " << fiberB << " = " << abs(fiber-fiberB) << std::endl;

	if( abs(fiber-fiberA)==1 ||
	    abs(fiber-fiberB)==1 ){
	  hitC=hit; ++flag[j]; break;
	}
	// if( (abs(fiber-fiberA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
	//     (abs(fiber-fiberB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
	//   hitC=hit; ++flag[j]; break;
	// }
      }
      if(hitC){
	TrCluster *cluster=new TrCluster(hitA,hitB,hitC);
	if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
	TrCluster *cluster=new TrCluster(hitA,hitB);
	if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      TrCluster *cluster=new TrCluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}			


bool TrAnalyzer::TrackSearchBFT( void )
{  
   //std::cout << "******************" << std::endl;
   
   // ConfMan *confMan = ConfMan::GetConfManager();
   // const int MinNumOfHitsTr = confMan->MinLayer();

   clearTracksBFT();
    
   int ntrack =
      LocalTrackSearch( &(BFTCHC[1]), TrackBFTCol, 
                        NumOfLayersBFT, MinNumOfHitsBFT );
   
  return true;
}

void TrAnalyzer::clearBFTHits( void )
{
   for( int l=0; l<=NumOfLayersBFT; ++l ){
      for_each( BFTHC[l].begin(),  BFTHC[l].end(),  DeleteObject() );
      BFTHC[l].clear();
   }
}

void TrAnalyzer::clearBFTCHits( void )
{
   for( int l=0; l<=NumOfLayersBFT; ++l ){
      for_each( BFTCHC[l].begin(),  BFTCHC[l].end(),  DeleteObject() );
      BFTCHC[l].clear();
      //hitlayersCont.clear();
   }
}

void TrAnalyzer::clearTracksBFT( void )
{
  for_each( TrackBFTCol.begin(), TrackBFTCol.end(), DeleteObject() );
  TrackBFTCol.clear();
}

bool TrAnalyzer::ReCalcBFTHits( bool applyRecursively )
{
  for( int l=0; l<=NumOfLayersBFT; ++l ){
    int n=BFTHC[l].size();
    for( int i=0; i<n; ++i ){
      TrHit *hit=(BFTHC[l])[i];
      if(hit) hit->ReCalcTr(applyRecursively);
    }
  }

  return true;
}

bool TrAnalyzer::ReCalcBFTClusters( bool applyRecursively )
{
  int n=BFTClCont.size();
  for( int i=0; i<n; ++i ){
    TrCluster *cl=BFTClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool TrAnalyzer::ReCalcTrackBFT( bool applyRecursively )
{
  int n=TrackBFTCol.size();
  for( int i=0; i<n; ++i ){
    TrLocalTrack *track=TrackBFTCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}



bool TrAnalyzer::ReCalcAll( void )
{
  ReCalcBFTHits();
  ReCalcTrackBFT();
  ReCalcBFTClusters();

  return true;
}

