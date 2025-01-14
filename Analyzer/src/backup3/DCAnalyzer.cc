/*
  DCAnalyzer.cc
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
#include "K18Track.hh"
#include "SksTrack.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "TemplateLib.hh"
#include "DCTrackSearch.hh"
#include "K18Parameters.hh"
#include "ConfMan.hh"
#include "MWPCCluster.hh"

// #include "Hodo2Hit.hh"
// #include "Hodo1Hit.hh"

//Add 2010/09/28///////////
// #include "SimuData.hh"
// #include "SimuDetectorId.hh"
////////////////////////////

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

//#include "SksObjectId.hh"

const double MaxChiSqrSksTrack = 10000.;
const double MaxTimeDifMWPC = 100.;

const double kMWPCClusteringWireExtension = 1.0;  // [mm]
const double kMWPCClusteringTimeExtension = 10.0; // [nsec]

const double Deg2Rad = acos(-1.)/180. ;

#ifdef MemoryLeak
debug::Counter DCAnalyzer::sm_counter("DCAnalyzer");
#endif

/******** Bc Out Tracking  *****/
//ver.2 -> by using VXU Tracking
#define ver2_BcOut  1

//ver.1 -> Nowadays Tracking
#define ver1_BcOut  0 
/*******************************/


/******** Sdc In Tracking  *****/
//ver.2 -> by using VXU Tracking
#define ver2_SdcIn  1 

//ver.1 -> Nowadays Tracking
#define ver1_SdcIn  0 
/*******************************/


  ///////////////////////////
  //SKS Tracking Rootine//
  ///////////////////////////
  
  //ver.2 -> New SksTracking Rootine
#define ver2_Sks  0

  //ver.1 -> Old SksTracking Rootine
#define ver1_Sks  1 

  
  double PIni = 0.90 ;//Initial Momentum (GeV/c)
  //  double MinP = 0.6 ;//300A Minimum Momentum
  double MinP = 0.7 ;//400A Minimum Momentum

  //measure the Initial Momentum by using Kinematics
#define Pini_kine  0 

  // Magnetic Field Strength of SKS
  double B_Sks = 2.1 ;//(T)
  
  ///////////////////////////
  ///////////////////////////

//______________________________________________________________________________
// forward declaration
int
clusterizeMWPCHit(const DCHitContainer&,
 		  MWPCClusterContainer& clusters);

//______________________________________________________________________________

DCAnalyzer::DCAnalyzer()
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

DCAnalyzer::~DCAnalyzer()
{
  clearSksTracks();
  clearK18Tracks();
  clearTracksSdcOut();
  clearTracksSdcIn();
  clearTracksBcOut();
  clearTracksBcIn();

  clearTracksBcOutSdcIn(); 

  clearDCHits();
  /*New Added by K.Miwa*/
  clearVtxHits();

#ifdef MemoryLeak
  --sm_counter;
#endif
}

bool DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[DCAnalyzer::DecodeRawHits]";

  clearDCHits();

  // BcIn
  for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *thit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	thit->SetTdcVal( rhit->GetTdc(j) );
	thit->SetTdcTrailing( rhit->GetTrailing(j) );
      }
      if(!thit) continue;
    
      if(thit->CalcMWPCObservables())
	TempBcInHC[layer].push_back(thit);
      else
	delete thit;
    }

    //std::cout<<"*************************************"<<std::endl;
    int ncl = clusterizeMWPCHit( TempBcInHC[layer], MWPCClCont[layer]);
    //std::cout<<"numCl="<< ncl << std::endl;
     for( int i=0; i<ncl; ++i ){
      MWPCCluster *p=MWPCClCont[layer][i];
      if (!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      double mwire    = mean.m_wire;
      double mwirepos = mean.m_wpos;
      double mtime    = mean.m_leading;
      double mtrail   = mean.m_trailing;
      
      DCHit *hit=new DCHit( layer+PlOffsBc, mwire );
      hit->SetClusterSize( p->GetClusterSize() );
      hit->SetMWPCFlag( true );
      hit->SetWire( mwire );
      hit->SetMeanWire( mwire );
      hit->SetMeanWirePosition( mwirepos );
      hit->SetDriftTime( mtime );
      hit->SetTrailingTime( mtrail );
      hit->SetDriftLength( 0. );
      hit->SetTdcVal( 0 );
      hit->SetTdcTrailing( 0 );
      if(!hit) continue;

      if(hit->CalcMWPCObservables())
	BcInHC[layer].push_back(hit);
      else
	delete hit;
    }
    //std::cout << "nh="<< BcInHC[layer].size() <<std::endl;
  } 
  
  // BcOut
  for( int layer=1; layer<=NumOfLayersBcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetBcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];
      
      DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
       int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;
      
      if(hit->CalcDCObservables())
	BcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  // SdcIn
  for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcInRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue; 

      if(hit->CalcDCObservables())
	SdcInHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  // SdcOut
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCRHitContainer &cont=rawData->GetSdcOutRawHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit=cont[i];

      DCHit *hit=new DCHit( rhit->PlaneId(), rhit->WireId() );
      int nhtdc= rhit->GetTdcSize();
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }
      if(!hit) continue;

      if(hit->CalcDCObservables())
	SdcOutHC[layer].push_back(hit);
      else
	delete hit;
    }
  }
      
  return true;
}

/* New Added by K.Miwa */
// bool DCAnalyzer::DecodeSimuHits( SimuData *simuData )
// {
//   const std::string funcname = "[DCAnalyzer::DecodeSimuHits]";

//   clearDCHits();
// #if 0
//   // BcIn
//   for( int layer=1; layer<=NumOfLayersBcIn; ++layer ){
//     const DCRHitContainer &cont=rawData->GetBcInRawHC(layer);
//     int nh=cont.size();
//     for( int i=0; i<nh; ++i ){
//        DCRawHit *rhit=cont[i];
//        if(!rhit) continue;
//        DCHit *hit=new DCHit(rhit->PlaneId()+PlOffsBc,
// 			    rhit->WireId(),rhit->Tdc());
//        if(!hit) continue;

//        if(hit->CalcObservables())
// 	 BcInHC[layer].push_back(hit);
//        else
// 	 delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
//  #endif

// #if 0
//   // BcOut
//   for( int layer=SimuPlMinBdcOut+SimuPlOffsBdc; layer<=SimuPlMaxBdcOut+SimuPlOffsBdc; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//       //      DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i), -1);
//       if(!hit) continue;
      
//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
//  	BcOutHC[layer-(SimuPlMinBdcOut+SimuPlOffsBdc-1)].push_back(hit);
//        else
// 	 delete hit;
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
// #endif  



//   // SdcIn
//   for( int layer=1; layer<=10; ++layer ){
//     int nh=simuData->GetDcMulti(layer);
//     for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//       //      DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1);
//       if(!hit) continue;
//       hit->SetTdcVal(-1);

//       if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i))){
//  	SdcInHC[layer].push_back(hit);
//       } 
//       else {
// 	delete hit;
//       }
//     } /* for( int i... ) */
//   } /* for( int layer... ) */
      

//   // SdcOut
//   for( int layer=31; layer<=42; ++layer ){
//      int nh=simuData->GetDcMulti(layer);
//      for( int i=0; i<nh; ++i ){
//       //Add 2010/09/29
//       DCHit *hit=new DCHit(layer, simuData->GetDcHitWire(layer, i));
//        //       DCHit *hit=new DCHit(layer,simuData->GetDcHitWire(layer, i),-1); 
//        if(!hit) continue;
//        hit->SetTdcVal(-1);
       
//        if(hit->CalcObservablesSimulation(simuData->GetDcDrLength(layer, i)))
// 	 SdcOutHC[layer-30].push_back(hit);
//        else
// 	 delete hit;
//      } /* for( int i... ) */
//   } /* for( int layer... ) */
  

//   return true;
// }

bool DCAnalyzer::TrackSearchBcIn( void )
{
  int ntrack =
    MWPCLocalTrackSearch(  &(BcInHC[1]), TrackBcInCol );

  return true;
}

bool DCAnalyzer::TrackSearchBcIn( const std::vector<std::vector<DCHitContainer> >& hc)
{
//   std::cout << __FILE__ << ":" << __func__ << std::endl;
  int ntracks = 0;
  MWPCLocalTrackSearch(hc, TrackBcInCol);

  std::vector<DCHitContainer> unusedHits;
  int toSearch=findUnusedHits(NumOfLayersBcIn, &(BcInHC[1]), unusedHits);
  if (toSearch>0) 
    {
//       std::cout << "#D DCAnalyzer::" << __func__ 
// 		<< " : " << ntracks << " tracks found. "
// 		<< " trying to find tracks w/o filter" << std::endl;
      std::vector<DCLocalTrack*> additionalTracks;
      int toAdd = MWPCLocalTrackSearch(&(unusedHits[0]), additionalTracks);
//       if (toAdd>0)
// 	std::cout << "#D " << __func__
// 		  << " : " << toAdd << " tracks found" << std::endl;
      TrackBcInCol.insert(TrackBcInCol.end(), 
			  additionalTracks.begin(), additionalTracks.end());
      //       ntracks += toAdd
    }
//   std::cout << "#D " << __FILE__ << ":" << __LINE__
// 	    << " " << __func__ << " total tracks = "
// 	    << TrackBcInCol.size() << std::endl;
  return true;
}


// Nowadays Tracking Routine
#if ver1_BcOut
bool DCAnalyzer::TrackSearchBcOut( void )
{
  int ntrack =
    LocalTrackSearch( BcOutHC, PPInfoBcOut, NPPInfoBcOut,
		      TrackBcOutCol, MinNumOfHitsBcOut );

  return true;
}
#endif

//   Add Y.Yonemoto 2010/6/25 (BcOut & SdcIn) by using VXU Tracking 
#if ver2_BcOut         
bool DCAnalyzer::TrackSearchBcOut( void )
{
  int ntrack =
    LocalTrackSearchVUX( BcOutHC, PPInfoBcOut, NPPInfoBcOut,
			 TrackBcOutCol, MinNumOfHitsBcOut );
  
  return true;
}
#endif


// Nowadays Tracking Routine
#if ver1_SdcIn
bool DCAnalyzer::TrackSearchSdcIn( void )
{
  int ntrack =
    LocalTrackSearch( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
		      TrackSdcInCol, MinNumOfHitsSdcIn );

  return true;
}
#endif

//   Add Y.Yonemoto 2010/6/25 (BcOut & SdcIn) by using VXU Tracking 
#if ver2_SdcIn         
bool DCAnalyzer::TrackSearchSdcIn( void )
{
  int ntrack =
    LocalTrackSearchVUX( SdcInHC, PPInfoSdcIn, NPPInfoSdcIn,
			 TrackSdcInCol, MinNumOfHitsSdcIn );
  
  //std::cout << "*************************" << std::endl;
  return true;
}
#endif

bool DCAnalyzer::TrackSearchSdcOut( void )
{

  //  int ntrack =
  //    SdcOutLocalTrackSearch( &(SdcOutHC[1]), TrackSdcOutCol, TrackSdcOutCol1, TrackSdcOutCol2);

  int ntrack =
    SdcOutLocalTrackSearch( SdcOutHC, TrackSdcOutCol );

  return true;
}

bool DCAnalyzer::TrackSearchBcOutSdcIn( void )
{

  int ntrack =
    LocalTrackSearchBcOutSdcIn( BcOutHC, PPInfoBcOut,
				SdcInHC, PPInfoSdcIn, 
				NPPInfoBcOut,NPPInfoSdcIn,
				TrackBcOutSdcInCol, 19 );

  return true;
}

bool DCAnalyzer::TrackSearchK18( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchK18]";

  clearK18Tracks();
  int nIn=TrackBcInCol.size(), nOut=TrackBcOutCol.size();
#if 0
  std::cout<<"**************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
	    << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif

  if( nIn==0 || nOut==0 ) return true;

  double pK18=ConfMan::GetConfManager()->K18Momentum();

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackBcInCol[iIn];
#if 0
    std::cout << "TrackIn  :" << std::setw(2) << iIn
	      << " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
	      << " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
	      << std::endl;
#endif
    if( !trIn->GoodForTracking() ||
	trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
	trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
	trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
	trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackBcOutCol[iOut];
#if 0
      std::cout << "TrackOut :" << std::setw(2) << iOut
		<< " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
		<< " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
		<< std::endl;
#endif
      if( !trOut->GoodForTracking() ||
	  trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
	  trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
	  trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
	  trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV ) continue;

#if 0
      std::cout << funcname << ": In -> " << trIn->GetChiSquare() 
		<< " (" << std::setw(2) << trIn->GetNHit() << ") "
		<< "Out -> " << trOut->GetChiSquare()
		<< " (" << std::setw(2) << trOut->GetNHit() << ") "
		<< std::endl;
#endif

      K18Track *track=new K18Track( trIn, trOut, pK18 );
      if( track && track->doFit() )
	K18TrackCol.push_back(track);
      else
	delete track;
    }
  }

#if 0
  std::cout<<"********************"<<std::endl;
 {
   int nn=K18TrackCol.size();
   std::cout << funcname << ": Before sorting. #Track=" 
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     K18Track *tp=K18TrackCol[i];

     std::cout << std::setw(3) << i 
	       << " ChiSqr=" << tp->chisquare()
	       << " Delta=" << tp->Delta()
	       << " P=" << tp->P() << "\n";
//      std::cout<<"********************"<<std::endl;
//      std::cout << "In :"
// 	       << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
// 	       << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
// 	       << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
// 	       << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
// 	       << "\n";
//      std::cout << "Out:"
// 	       << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
// 	       << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
// 	       << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
// 	       << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
// 	       << std::endl;
   }
 }
#endif

  partial_sort( K18TrackCol.begin(), K18TrackCol.end(),
		K18TrackCol.end(), K18TrackComp() );

  return true;
}


bool DCAnalyzer::TrackSearchSks( void )
{
  static const std::string funcname = "[DCAnalyzer::TrackSearchSks]";

  clearSksTracks();

  int nIn=TrackSdcInCol.size(), nOut=TrackSdcOutCol.size();
#if 0
  std::cout<<"*********************************************"<<std::endl;
  std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn
	    << " #TracksOut=" << std::setw(3) << nOut << std::endl;
#endif
  if( nIn==0 || nOut==0 ) return true;

  for( int iIn=0; iIn<nIn; ++iIn ){
    DCLocalTrack *trIn=TrackSdcInCol[iIn];
    if( !trIn->GoodForTracking() ) continue;
    for( int iOut=0; iOut<nOut; ++iOut ){
      DCLocalTrack *trOut=TrackSdcOutCol[iOut];
      if( !trOut->GoodForTracking() ) continue;


#if Pini_kine  
      /*************************************************************************/
      /*****   measure the Initial Momentum of SKS by using Kinematics   *******/

      double X0In = trIn->GetX0(), Y0In = trIn->GetY0();
      double U0In = trIn->GetU0(), V0In = trIn->GetV0();
      double X0Out = trOut->GetX0(), Y0Out = trOut->GetY0();
      double U0Out = trOut->GetU0(), V0Out = trOut->GetV0();
      
      int SksVp1_Id = DCGeomMan::GetInstance().GetSksVp1Id();
      int SksVp5_Id = DCGeomMan::GetInstance().GetSksVp5Id();
      double zSksVp1 = DCGeomMan::GetInstance().GetLocalZ( SksVp1_Id );
      double zSksVp5 = DCGeomMan::GetInstance().GetLocalZ( SksVp5_Id );
      
  
      //Get Global Position at the SksVp1 & SksVp5
      ThreeVector gpVp1 = DCGeomMan::GetInstance().GetGlobalPosition(SksVp1_Id);
      ThreeVector gpVp5 = DCGeomMan::GetInstance().GetGlobalPosition(SksVp5_Id);
      
      
      double X0_Vp1 = gpVp1.x(); 
      double Y0_Vp1 = gpVp1.y(); 
      double X0_Vp5 = gpVp5.x(); 
      double Y0_Vp5 = gpVp5.y(); 
      
      
      //SdcIn Local Tracking -> X(Z = SksVp1) = X0 + U0*Z 
      double XIn_Vp1 = X0In + U0In*zSksVp1 ;
      double theta_In = atan(U0In);
      
      //triangle 
      double dX_Vp1 = fabs(XIn_Vp1);
      
      double L_Vp1 = dX_Vp1*sin((90-theta_In)*Deg2Rad)/sin((60+theta_In)*Deg2Rad);
      
      double X_Vp1,Y_Vp1;
      
      if(XIn_Vp1>=0){
	X_Vp1 = X0_Vp1 - L_Vp1*cos(70*Deg2Rad);
	Y_Vp1 = Y0_Vp1 - L_Vp1*sin(70*Deg2Rad);
      }
      if(XIn_Vp1<0){
	X_Vp1 = X0_Vp1 + L_Vp1*cos(70*Deg2Rad);
	Y_Vp1 = Y0_Vp1 + L_Vp1*sin(70*Deg2Rad);
      }
      
      //SdcOut Local Tracking -> X(Z = SksVp5) = X0 + U0*Z 
      double XOut_Vp5 = X0Out + U0Out*zSksVp5 ;
      double theta_Out = atan(U0Out);
      
      //triangle 
      double dX_Vp5 = fabs(XOut_Vp5);
      double L_Vp5 = dX_Vp5*sin((90+theta_Out)*Deg2Rad)/sin((60-theta_Out)*Deg2Rad);
      
      double X_Vp5,Y_Vp5;
      
      if(XOut_Vp5>=0){
	X_Vp5 = X0_Vp5 + L_Vp5*cos(70*Deg2Rad);
	Y_Vp5 = Y0_Vp5 - L_Vp5*sin(70*Deg2Rad);
      }
      if(XOut_Vp5<0){
	X_Vp5 = X0_Vp5 - L_Vp5*cos(70*Deg2Rad);
	Y_Vp5 = Y0_Vp5 + L_Vp5*sin(70*Deg2Rad);
      }
      
      
      double L_Sks = sqrt(pow( (X_Vp1-X_Vp5) ,2) 
			  + pow( (Y_Vp1-Y_Vp5),2))  ;
      
      double theta_Sks = 100 - theta_In + theta_Out ;
      
      double r_Sks = L_Sks/2/sin(theta_Sks*Deg2Rad/2);
      
      double P_xy = 0.3 * B_Sks * r_Sks *0.001;
      
      //Momentum of SKS by using Kinematics 
      double P_Sks = P_xy*sqrt(1+V0In*V0In+U0In*U0In)/sqrt(1+U0In*U0In) ;
      
      //Momentum Correction (Please Check the figure of the relationship of P_Sks & U0In) 
      double P_Sks_Cor = P_Sks+0.3204898*U0In + 0.12 ;
      
      //Initial momentum
      Pini = P_Sks_Cor ;

      /*****   measure the Initial Momentum of SKS by using Kinematics   *******/
      /*************************************************************************/
#endif



#if ver2_Sks

      SksTrack *pretrack = new SksTrack( trIn, trOut );      
      int preK=0,prekk=2;
      //double PIni = 1.1 ;
      //      double PIni = 0.85 ;
      //      double MinP = 0.85 ;
      double preiP,preiP_pre=PIni;
      int preNiteration=0;
      bool prestatus;
      double chisqr2,prep_check,prechisqr2=10000,prechisqr3=1.0e+10,preP2=PIni;
      double dchisqr2=-1;

      //      std::cout << "pre SKSTracking start!!" << std::endl;      

      for(;;){
	if(preK==0) preiP = PIni;
	pretrack->SetInitialMomentum( preiP );
	pretrack->PreTracking( 1 );//2010-7-28
	pretrack->doFit();
	preNiteration++;
	//	preNiteration=pretrack->Niteration();
	
	prep_check=pretrack->PrimaryMomentum().mag();
	prestatus=pretrack->Status();
	chisqr2=pretrack->chisqr();
	
	if((prechisqr2>chisqr2) && (prestatus)){
	  prechisqr2 = chisqr2 ;
	  preP2 = preiP ;
	}
	
#if 0
	  /*
	    std::cout << "Event # = " << Number << " : prestatus = " << prestatus << std::endl;
	    std::cout << "preK = " << preK << " : # of iteration = " << pretrack->Niteration() << " : preiP = " << preiP << " : P_g4 = " << gp << std::endl;
	    std::cout << "P = " << prep_check << " : chisqr2 = " << chisqr2  << std::endl;
	    std::cout << "preP2 = " << preP2 << " : prechisqr2 = " << prechisqr2  << std::endl;
	  */
	  //	  std::cout << pretrack->Niteration() << " " << chisqr2 << " " << preiP  << " " << PIni << std::endl;
#endif
	
	++preK;
	if((chisqr2>1.0e+05) || !(chisqr2>1.0e-100))
	  preiP=preiP-0.02;
	else if(chisqr2>1.0e+04)
	  preiP=preiP-0.01;
	else if(chisqr2>1.0e+03)
	  preiP=preiP-0.002;//0.001
	else if(chisqr2>1.0e+02)
	  preiP=preiP-0.001;//0.0005
	else 
	  preiP=preiP-0.0005;//0.0001

	prechisqr3 = chisqr2;

	//delete pretrack;
	if( ((log10(chisqr2)-log10(prechisqr2))>1) && (prechisqr2<1.0e+02) && (prechisqr2>0) ) break;
	if(preiP<MinP) break;
      }
     
      if((prechisqr2<1.0e+03) && (prechisqr2>0)){
	/***** Sks Traking Start!! ************/
	SksTrack *track = new SksTrack( trIn, trOut );
	int K=0,kk=2;
	double iP,iP_pre=preP2;
	//double iP,iP_pre=0.75;
	//  double iP,iP_pre=preP2+0.02;
	int Niteration;
	bool status;
	double chisqr,p_check,prechisqr,preP;
	
	for(;;){
	  if(K==0) iP = iP_pre;
	  track->SetInitialMomentum( iP );
	  //    track->SetInitialMomentum( PIni );
	  track->doFit();
	  Niteration=track->Niteration();
	  p_check=track->PrimaryMomentum().mag();
	  status=track->Status();
	  chisqr=track->chisqr();
	  
	  if(K==0){
	    prechisqr = chisqr ;
	    preP = iP ;
	  }
	  if(K>0){
	    if(prechisqr>chisqr){
	      prechisqr = chisqr ;
	      preP = iP ;
	    }
	  }
	  
#if 1
	  std::cout << K << " " << track->Niteration() << " " << chisqr 
		    << " " << p_check << " " << preP2 << std::endl;
#endif
	  
	  if(status || (p_check>0 && p_check<=MinP) || (iP<=MinP)) break;
	  if((K==0) && !(chisqr<0) && (p_check>MinP)){
	    iP_pre=p_check-0.01;
	  }
	  if(K!=0 && !(chisqr<0) && (p_check>MinP)){
	    iP=p_check-0.01;
	    if((iP_pre-iP)>0.01){
	      iP_pre=p_check-0.01;
	      kk = 2;
	    }
	    if(((iP_pre-iP)<=0.01)){ 
	      if(iP>0) iP=p_check-0.01*kk ;
	      else iP=iP_pre-0.01*kk ;
	      ++kk;
	    }
	  }
	  if(chisqr<0 || !(p_check>0)){
	    iP_pre=iP-0.01;
	    iP=iP-0.01;
	  }
	  
	  ++K;
	  //  iP=iP-0.005;
	  if(iP<MinP) break;
	}
	
	//std::cout << "chisqr = " << chisqr << " : momentum = " << p_check << std::endl;
	
	//std::cout << "SKSTracking finish!! ( # iteration = " << K << " )" << std::endl;      

	if(!track) continue;
	if( track->doFit() && track->chisqr()<MaxChiSqrSksTrack ){
	  SksTrackCol.push_back(track);
	}
	else{
	  delete track;
	}
      }
      delete pretrack;

#endif
      



#if ver1_Sks
      SksTrack *tp = new SksTrack( trIn, trOut );

      if(!tp) continue;
      if( tp->doFit() && tp->chisqr()<MaxChiSqrSksTrack ){
	SksTrackCol.push_back(tp);
      }
      else{
	delete tp;
      }
#endif
      
    }
    
  }






  partial_sort( SksTrackCol.begin(), SksTrackCol.end(),
		SksTrackCol.end(), SksTrackComp() );


#if 0
  std::cout<<"********************"<<std::endl;
 {
   int nn=SksTrackCol.size();
   std::cout << funcname << ": Before Deleting. #Track="
	     << nn << std::endl;
   for( int i=0; i<nn; ++i ){
     SksTrack *tp=SksTrackCol[i];
     std::cout << std::setw(3) << i
	       << " Nitra=" << std::setw(3) << tp->Niteration()
	       << " ChiSqr=" << tp->chisqr()
	       << " P=" << tp->PrimaryMomentum().mag()
	       << " PL(TOF)=" << tp->PathLengthToTOF()
	       << std::endl;
   }
 }
#endif

  return true;
}

void DCAnalyzer::clearDCHits( void )
{
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    for_each( TempBcInHC[l].begin(),  TempBcInHC[l].end(),  DeleteObject() );
    TempBcInHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    for_each( BcInHC[l].begin(),  BcInHC[l].end(),  DeleteObject() );
    BcInHC[l].clear();
  }

  for ( int l=0; l<=NumOfLayersBcIn; ++l ){
    std::for_each( MWPCClCont[l].begin(),  MWPCClCont[l].end(), DeleteObject() );
    MWPCClCont[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcOut; ++l ){
    for_each( BcOutHC[l].begin(), BcOutHC[l].end(), DeleteObject() );
    BcOutHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcIn; ++l ){
    for_each( SdcInHC[l].begin(),  SdcInHC[l].end(),  DeleteObject() );
    SdcInHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcOut; ++l ){
    for_each( SdcOutHC[l].begin(), SdcOutHC[l].end(), DeleteObject() );
    SdcOutHC[l].clear();
  }
}

void DCAnalyzer::clearVtxHits( void )
{
  for_each( VtxPoint.begin(),  VtxPoint.end(),  DeleteObject() );
  VtxPoint.clear();
}

void DCAnalyzer::clearTracksBcIn( void )
{
  for_each( TrackBcInCol.begin(), TrackBcInCol.end(), DeleteObject() );
  TrackBcInCol.clear();
}

void DCAnalyzer::clearTracksBcOut( void )
{
  for_each( TrackBcOutCol.begin(), TrackBcOutCol.end(), DeleteObject() );
  TrackBcOutCol.clear();
}

void DCAnalyzer::clearTracksSdcIn( void )
{
  for_each( TrackSdcInCol.begin(), TrackSdcInCol.end(), DeleteObject() );
  TrackSdcInCol.clear();
}

void DCAnalyzer::clearTracksSdcOut( void )
{
  for_each( TrackSdcOutCol.begin(), TrackSdcOutCol.end(), DeleteObject() );
  TrackSdcOutCol.clear();
}

void DCAnalyzer::clearK18Tracks( void )
{
  for_each( K18TrackCol.begin(), K18TrackCol.end(), DeleteObject() );
  K18TrackCol.clear();
}

void DCAnalyzer::clearSksTracks( void )
{
  for_each( SksTrackCol.begin(), SksTrackCol.end(), DeleteObject() );
  SksTrackCol.clear();
}

void DCAnalyzer::clearTracksBcOutSdcIn( void )
{
  for_each( TrackBcOutSdcInCol.begin(), TrackBcOutSdcInCol.end(),
            DeleteObject() );
  TrackBcOutSdcInCol.clear();
}

bool DCAnalyzer::ReCalcDCHits( bool applyRecursively )
{
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    int n=TempBcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(TempBcInHC[l])[i];
      if(hit) hit->ReCalcMWPC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersBcIn; ++l ){
    int n=BcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(BcInHC[l])[i];
      if(hit) hit->ReCalcMWPC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersBcOut; ++l ){
    int n=BcOutHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(BcOutHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersSdcIn; ++l ){
    int n=SdcInHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(SdcInHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }
  for( int l=0; l<=NumOfLayersSdcOut; ++l ){
    int n=SdcOutHC[l].size();
    for( int i=0; i<n; ++i ){
      DCHit *hit=(SdcOutHC[l])[i];
      if(hit) hit->ReCalcDC(applyRecursively);
    }
  }

  return true;
}

bool DCAnalyzer::ReCalcTrackBcIn( bool applyRecursively )
{
  int n=TrackBcInCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackBcInCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackBcOut( bool applyRecursively )
{
  int n=TrackBcOutCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackBcOutCol[i];
    if( track ) track->ReCalc( applyRecursively ); 
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackSdcIn( bool applyRecursively )
{
  int n=TrackSdcInCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackSdcInCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcTrackSdcOut( bool applyRecursively )
{
  int n=TrackSdcOutCol.size();
  for( int i=0; i<n; ++i ){
    DCLocalTrack *track=TrackSdcOutCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcK18Track( bool applyRecursively )
{
  int n=K18TrackCol.size();
  for( int i=0; i<n; ++i ){
    K18Track *track=K18TrackCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcSksTrack( bool applyRecursively )
{
  int n=SksTrackCol.size();
  for( int i=0; i<n; ++i ){
    SksTrack *track=SksTrackCol[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

bool DCAnalyzer::ReCalcAll( void )
{
  ReCalcDCHits();

  ReCalcTrackBcIn();
  ReCalcTrackBcOut();
  ReCalcTrackSdcIn();
  ReCalcTrackSdcOut();

  ReCalcK18Track();
  ReCalcSksTrack();

  return true;
}

//______________________________________________________________________________
double
log2(double x)
{
  return std::log(x)/std::log(2.0);
}

//______________________________________________________________________________
bool
isConnectable(double wire1,
	      double leading1,
	      double trailing1,
	      double wire2,
	      double leading2,
	      double trailing2,
	      double wExt,
	      double tExt)
{
  double w1Min = wire1 - wExt;
  double w1Max = wire1 + wExt;
  double t1Min = leading1  - tExt;
  double t1Max = trailing1 + tExt;

  double w2Min = wire2 - wExt;
  double w2Max = wire2 + wExt;
  double t2Min = leading2  - tExt;
  double t2Max = trailing2 + tExt;

  bool isWireOk = !(w1Min>w2Max || w1Max<w2Min);
  bool isTimeOk = !(t1Min>t2Max || t1Max<t2Min);

//   std::cout << " w1 = " << wire1
// 	    << " le1 = " << leading1
// 	    << " tr1 = " << trailing1 << "\n"
// 	    << " w2 = " << wire2
// 	    << " le2 = " << leading2
// 	    << " tr2 = " << trailing2 << "\n"
// 	    << " w1(" << w1Min << " -- " << w1Max << "), t1(" 
// 	    << t1Min << " -- " << t1Max << ")\n" 
// 	    << " w2(" << w2Min << " -- " << w2Max << "), t2(" 
// 	    << t2Min << " -- " << t2Max << ")\n" 
// 	    << " wire : " << isWireOk
// 	    << ", time : " << isTimeOk
// 	    << std::endl;

  return (isWireOk && isTimeOk);
}

//______________________________________________________________________________
void
printConnectionFlag(const std::vector<std::deque<bool> >& flag)
{
  for (int i=0, n=flag.size(); i<n; ++i)
    {
      std::cout << "\n";
      for (int j=0, m=flag[i].size(); j<m; ++j)
	{
	  std::cout << " " << flag[i][j];
	}
    }
  std::cout << std::endl;
  return;
}

//______________________________________________________________________________
int
clusterizeMWPCHit(const DCHitContainer& hits,
		  MWPCClusterContainer& clusters)
{
  if (!clusters.empty())
    {
      std::for_each(clusters.begin(), clusters.end(), DeleteObject());
      clusters.clear();
    }


  const int nhits = hits.size();
//   std::cout << "#D " << __func__ << " " << nhits << std::endl;
  if (nhits==0)
    return 0;

  int n = 0;
  for (int i=0; i<nhits; ++i)
    {
      const DCHit* h = hits[i];
      if (!h)
	continue;
      n += h->GetTdcSize();
    }

  DCHitContainer singleHits;
  singleHits.reserve(n);
  for (int i=0; i<nhits; ++i)
    {
      const DCHit* h = hits[i];
      if (!h)
	continue;
      int nn = h->GetTdcSize();
      for (int ii=0; ii<nn; ++ii)
	{
	  DCHit* htmp = new DCHit(h->GetLayer(), h->GetWire());
	  htmp->SetTdcVal(h->GetTdcVal());
	  htmp->SetTdcTrailing(h->GetTdcTrailing());
	  htmp->SetTrailingTime(h->GetTrailingTime());
	  htmp->SetDriftTime(h->GetDriftTime());
	  htmp->SetDriftLength(h->GetDriftLength());
	  htmp->SetTiltAngle(h->GetTiltAngle());
	  htmp->SetWirePosition(h->GetWirePosition());
	  htmp->setRangeCheckStatus(h->rangecheck(), 0);
	  singleHits.push_back(htmp);
	}
    }
  
  std::vector<std::deque<bool> > flag(n, std::deque<bool>(n, false));
  n = singleHits.size();
  for (int i=0;  i<n; ++i)
    {
      flag[i][i] = true;
      const DCHit* h1 = singleHits[i];
//       h1->print("h1");
      for (int j=i+1; j<n; ++j)
	{
	  const DCHit* h2 = singleHits[j];
// 	  h2->print("h2");
// 	  std::cout << " (i,j) = (" << i << ", " << j << ")" << std::endl;
	  bool val
	    = isConnectable(h1->GetWirePosition(), 
			    h1->GetDriftTime(),
			    h1->GetTrailingTime(),
			    h2->GetWirePosition(),
			    h2->GetDriftTime(),
			    h2->GetTrailingTime(),
			    kMWPCClusteringWireExtension,
			    kMWPCClusteringTimeExtension);
// 	  std::cout << "#D val = " << val << std::endl;
	  flag[i][j] = val;
	  flag[j][i] = val;
	}
    }

//   std::cout << "#D " << __func__ << "  before " << std::endl;
//   printConnectionFlag(flag);

  const int maxLoop = static_cast<int>(log2(n))+1;
  for (int loop=0; loop<maxLoop; ++loop)
    {
      std::vector<std::deque<bool> > tmp(n, std::deque<bool>(n, false));
      for (int i=0; i<n; ++i)
	{
	  for (int j=i; j<n; ++j)
	    {
	      for (int k=0; k<n; ++k)
		{
		  tmp[i][j] |= (flag[i][k] && flag[k][j]);
		  tmp[j][i] = tmp[i][j];
		}
	    }
	}
      flag = tmp;
//       std::cout << " n iteration = " << loop << std::endl;
//       printConnectionFlag(flag);
    }

//   std::cout << "#D " << __func__ << "  after " << std::endl;
//   printConnectionFlag(flag);

  std::set<int> checked;
  for (int i=0; i<n; ++i)
    {
      if (checked.find(i)!=checked.end())
	continue;
      MWPCCluster* c = 0;
      for (int j=i; j<n; ++j)
	{
	  if (flag[i][j])
	    {
	      checked.insert(j);
	      if (!c) 
		{
		  c = new MWPCCluster;
// 		  std::cout << " new cluster " << std::endl;
		}
// 	      std::cout << " " << i << "---" << j << std::endl;
	      c->Add(singleHits[j]);
	    }
	}
      
      if (c)
	{
	  c->Calculate();
	  clusters.push_back(c);
	}
    }

//   std::cout << " end of " << __func__
// 	    << " : n = " << n << ", " << checked.size()
// 	    << std::endl;

//   std::cout << __func__ << " n clusters = " << clusters.size() << std::endl;
  return clusters.size();
}


////////////
// for DS //
////////////

// std::size_t DCAnalyzer::DSSize( void ) const
// {
//   static const std::string funcname = "[DCAnalyzer::DSSize]";
//   std::size_t totalSize=1;

// #if 0
//   std::cout << funcname << ": " << totalSize;
// #endif

//   // DCHit
//   for( int l=1; l<=NumOfLayersBcIn; ++l ){
//     int nh=BcInHC[l].size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += (BcInHC[l])[i]->DSSize();
//   }
//   for( int l=1; l<=NumOfLayersBcOut; ++l ){
//     int nh=BcOutHC[l].size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += (BcOutHC[l])[i]->DSSize();
//   }
//   for( int l=1; l<=NumOfLayersSdcIn; ++l ){
//     int nh=SdcInHC[l].size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += (SdcInHC[l])[i]->DSSize();
//   }
//   for( int l=1; l<=NumOfLayersSdcOut; ++l ){
//     int nh=SdcOutHC[l].size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += (SdcOutHC[l])[i]->DSSize();
//   }

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   // DCLocalTrack
//   {
//     int nh=TrackBcInCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += TrackBcInCol[i]->DSSize();
//   }
//   {
//     int nh=TrackBcOutCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += TrackBcOutCol[i]->DSSize();
//   }
//   {
//     int nh=TrackSdcInCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += TrackSdcInCol[i]->DSSize();
//   }
//   {
//     int nh=TrackSdcOutCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += TrackSdcOutCol[i]->DSSize();
//   }

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   // K18Track
//   {
//     int nh=K18TrackCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += K18TrackCol[i]->DSSize();
//   }

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   // SksTrack
//   {
//     int nh=SksTrackCol.size();
//     totalSize += 1;
//     for( int i=0; i<nh; ++i )
//       totalSize += SksTrackCol[i]->DSSize();
//   }

// #if 0
//   std::cout << " " << totalSize << std::endl; 
// #endif

//   return totalSize;
// }

// std::size_t DCAnalyzer::DSSave( unsigned int *bufp )
// {
//   static const std::string funcname = "[DCAnalyzer::DSSave]";
//   unsigned int *hp = bufp;
//   int HeaderID = DSMakeHeaderID( SksObjDCAna, true );
//   std::size_t totalSize = 1;
//   bufp += 1;
//   std::size_t size, size1;
//   unsigned int *hp1;
//   int nh;

// #if 0
//   std::cout << funcname << ": " << totalSize;
// #endif

//   ///////////
//   // DCHit //
//   ///////////
//   // BcIn
//   for( int l=1; l<=NumOfLayersBcIn; ++l ){
//     hp1=bufp; bufp+=1; size1=1; nh=BcInHC[l].size();
//     int header=DSMakeHeaderID( SksObjBcInHC, true, true, l );
//     for( int i=0; i<nh; ++i ){
//       size=(BcInHC[l])[i]->DSSave(bufp); bufp+=size; size1+=size;
//       (BcInHC[l])[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }
//   // BcOut
//   for( int l=1; l<=NumOfLayersBcOut; ++l ){
//     hp1=bufp; bufp+=1; size1=1; nh=BcOutHC[l].size();
//     int header=DSMakeHeaderID( SksObjBcOutHC, true, true, l );
//     for( int i=0; i<nh; ++i ){
//       size=(BcOutHC[l])[i]->DSSave(bufp); bufp+=size; size1+=size;
//       (BcOutHC[l])[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }
//   // SdcIn
//   for( int l=1; l<=NumOfLayersSdcIn; ++l ){
//     hp1=bufp; bufp+=1; size1=1; nh=SdcInHC[l].size();
//     int header=DSMakeHeaderID( SksObjSdcInHC, true, true, l );
//     for( int i=0; i<nh; ++i ){
//       size=(SdcInHC[l])[i]->DSSave(bufp); bufp+=size; size1+=size;
//       (SdcInHC[l])[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }
//   // SdcOut
//   for( int l=1; l<=NumOfLayersSdcOut; ++l ){
//     hp1=bufp; bufp+=1; size1=1; nh=SdcOutHC[l].size();
//     int header=DSMakeHeaderID( SksObjSdcOutHC, true, true, l );
//     for( int i=0; i<nh; ++i ){
//       size=(SdcOutHC[l])[i]->DSSave(bufp); bufp+=size; size1+=size;
//       (SdcOutHC[l])[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   //////////////////
//   // DCLocalTrack //
//   //////////////////
//   // BcIn
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=TrackBcInCol.size();
//     int header=DSMakeHeaderID( SksObjLocalTrackBIC, true );
//     for( int i=0; i<nh; ++i ){
//       size=TrackBcInCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       TrackBcInCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }    
//   // BcOut
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=TrackBcOutCol.size();
//     int header=DSMakeHeaderID( SksObjLocalTrackBOC, true );
//     for( int i=0; i<nh; ++i ){
//       size=TrackBcOutCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       TrackBcOutCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }    
//   // SdcIn
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=TrackSdcInCol.size();
//     int header=DSMakeHeaderID( SksObjLocalTrackSIC, true );
//     for( int i=0; i<nh; ++i ){
//       size=TrackSdcInCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       TrackSdcInCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }    
//   // SdcOut
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=TrackSdcOutCol.size();
//     int header=DSMakeHeaderID( SksObjLocalTrackSOC, true );
//     for( int i=0; i<nh; ++i ){
//       size=TrackSdcOutCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       TrackSdcOutCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }    

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   /////////////
//   // K18Track //
//   /////////////
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=K18TrackCol.size();
//     int header=DSMakeHeaderID( SksObjK18TrackC, true );
//     for( int i=0; i<nh; ++i ){
//       size=K18TrackCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       K18TrackCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }

// #if 0
//   std::cout << " " << totalSize;
// #endif

//   //////////////
//   // SksTrack //
//   //////////////
//   {
//     hp1=bufp; bufp+=1; size1=1; nh=SksTrackCol.size();
//     int header=DSMakeHeaderID( SksObjSksTrackC, true );
//     for( int i=0; i<nh; ++i ){
//       size=SksTrackCol[i]->DSSave(bufp); bufp+=size; size1+=size;
//       SksTrackCol[i]->DSSetSaveKey( DSMakeKey( header, i ) );
//     }
//     *hp1 = DSMakeHeader( header, size1 );
//     totalSize += size1;
//   }

//   *hp = DSMakeHeader( HeaderID, totalSize );

// #if 0
//   std::cout << " " << totalSize << std::endl;
// #endif

//   return totalSize;
// } 

// bool DCAnalyzer::DSRestore( unsigned int *bufp )
// {
//   static const std::string funcname =
//     "[DCAnalyzer::DSRestore]";

//   clearDCHits();
//   clearTracksBcIn(); clearTracksBcOut();
//   clearTracksSdcIn(); clearTracksSdcOut();
//   clearK18Tracks(); clearSksTracks();
  
//   unsigned int *hp=bufp;
//   std::size_t totalSize = DSGetSize( *bufp );
// #if 0
//   std::cout << funcname << ": 0 Size=" << std::setw(4) << totalSize
// 	    << " ObjId=" << std::setw(3)
// 	    << DSGetObjectID(DSGetHeader(*hp)) << std::endl;
// #endif

//   bufp+=1;
//   while( bufp < hp+totalSize ){
//     unsigned int *hp1=bufp;
//     std::size_t size = DSGetSize( *hp1 );
//     unsigned int head = DSGetHeaderID( *hp1 );
//     int objId = DSGetObjectID( head );
// #if 0
//     std::cout << funcname << ": 1 Size=" << std::setw(4) << size
// 	      << " ObjId=" << std::setw(3) << objId << std::endl;
// #endif
//     int subId;
//     std::size_t size1;
//     bufp+=1;

//     switch( objId ){
//     case SksObjBcInHC:
//       subId=DSGetObjectSubID( head );
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjBcInHC SubId=" << subId
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	if(subId>0 && subId<=NumOfLayersBcIn ){
// 	  try{ 
// 	    DCHit *hit=new DCHit( bufp );
// 	    if( hit ) BcInHC[subId].push_back( hit );
// 	  }
// 	  catch( std::exception &e ){
// 	    std::cerr << funcname << ":SksObjBcInHC: " << e.what()
// 		      << " subId=" << subId << " Size=" << size1 
// 		      << " ContSize=" << BcInHC[subId].size() << std::endl;
// 	    throw;
// 	  }
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjBcOutHC:
//       subId=DSGetObjectSubID( head );
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjBcOutHC SubId=" << subId
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	if(subId>0 && subId<=NumOfLayersBcOut ){
// 	  try{ 
// 	    DCHit *hit=new DCHit( bufp );
// 	    if( hit ) BcOutHC[subId].push_back( hit );
// 	  }
// 	  catch( std::exception &e ){
// 	    std::cerr << funcname << ":SksObjBcOutHC: " << e.what()
// 		      << " subId=" << subId << " Size=" << size1 
// 		      << " ContSize=" << BcOutHC[subId].size() << std::endl;
// 	    throw;
// 	  }

// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjSdcInHC:
//       subId=DSGetObjectSubID( head );
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjSdcInHC SubId=" << subId
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	if(subId>0 && subId<=NumOfLayersSdcIn ){
// 	  try{
// 	    DCHit *hit=new DCHit( bufp );
// 	    if( hit ) SdcInHC[subId].push_back( hit );
// 	  }
// 	  catch( std::exception &e ){
// 	    std::cerr << funcname << ":SksObjSdcInHC: " << e.what()
// 		      << " subId=" << subId << " Size=" << size1 
// 		      << " ContSize=" << SdcInHC[subId].size() << std::endl;
// 	    throw;
// 	  }
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjSdcOutHC:
//       subId=DSGetObjectSubID( head );
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjSdcOutHC SubId=" << subId
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	if(subId>0 && subId<=NumOfLayersSdcOut ){
// 	  try{ 
// 	    DCHit *hit=new DCHit( bufp );
// 	    if( hit ) SdcOutHC[subId].push_back( hit );
// 	  }
// 	  catch( std::exception &e ){
// 	    std::cerr << funcname << ":SksObjSdcOutHC: " << e.what()
// 		      << " subId=" << subId << " Size=" << size1 
// 		      << " ContSize=" << SdcOutHC[subId].size() << std::endl;
// 	    throw;
// 	  }
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjLocalTrackBIC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjLocalTrackBIC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{ 
// 	  DCLocalTrack *tp=new DCLocalTrack( bufp, this );
// 	  if( tp ) TrackBcInCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjLocalTrackBIC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << TrackBcInCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjLocalTrackBOC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjLocalTrackBOC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{
// 	  DCLocalTrack *tp=new DCLocalTrack( bufp, this );
// 	  if( tp ) TrackBcOutCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjLocalTrackBOC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << TrackBcOutCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjLocalTrackSIC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjLocalTrackSIC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{
// 	  DCLocalTrack *tp=new DCLocalTrack( bufp, this );
// 	  if( tp ) TrackSdcInCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjLocalTrackSIC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << TrackSdcInCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjLocalTrackSOC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjLocalTrackSOC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{
// 	  DCLocalTrack *tp=new DCLocalTrack( bufp, this );
// 	  if( tp ) TrackSdcOutCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjLocalTrackSOC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << TrackSdcOutCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjK18TrackC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjK18TrackC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{ 
// 	  K18Track *tp=new K18Track( bufp, this );
// 	  if( tp ) K18TrackCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjK18TrackC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << K18TrackCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     case SksObjSksTrackC:
//       while( bufp<hp1+size ){
// 	size1=DSGetSize(*bufp);
// #if 0
// 	std::cout << funcname << ": SksObjSksTrackC " 
// 		  << "  Obj=" << DSGetObjectID( DSGetHeaderID(*bufp) ) 
// 		  << std::endl;
// #endif
// 	try{
// 	  SksTrack *tp=new SksTrack( bufp, this );
// 	  if( tp ) SksTrackCol.push_back(tp);
// 	}
// 	catch( std::exception &e ){
// 	  std::cerr << funcname << ":SksObjSkslTrackC: " << e.what()
// 		    << " Size=" << size1 << " ContSize=" << SksTrackCol.size() 
// 		    << std::endl;
// 	  throw;
// 	}
// 	bufp+=size1;
//       }
//       break;
//     default:
//       {
// 	std::ostringstream mess;
// 	mess << funcname << ": invalid object type ( " << objId  << " )";
// 	std::cerr << mess << std::endl;
// 	throw std::invalid_argument(mess.str());
//       }
//       break;
//     } /* swith( objId ) */
//   } /* while( ... ) */

//   return true;
// }
