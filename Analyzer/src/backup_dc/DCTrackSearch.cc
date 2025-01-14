/*
  DCTrackSearch.cc
*/

#include "DCTrackSearch.hh"
#include "DCParameters.hh"
#include "DCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"
#include "MWPCCluster.hh"

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double MaxChisquare = 10000.;
const double MaxNumberOfClusters = 30.;
const double MaxCombi = 1000000.;//1,000,000

//For BC3&4, SDC1&2
int LocalTrackSearch( const DCHitContainer * HC,
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> &TrackCont,
		      int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(npp);

  for( int i=0; i<npp; ++i ){
    bool ppFlag=PpInfo[i].flag;
    int layer1=PpInfo[i].id1, layer2=PpInfo[i].id2;
     if(ppFlag) 
      MakePairPlaneHitCluster( HC[layer1], HC[layer2], 
			       PpInfo[i].CellSize, CandCont[i] );
    else
      MakeUnPairPlaneHitCluster( HC[layer1], CandCont[i] );
  }

  std::vector <int> nCombi(npp);
  for( int i=0; i<npp; ++i ){ 
    nCombi[i]=(CandCont[i]).size();
    if( nCombi[i]>MaxNumberOfClusters ) nCombi[i]=0;
  }
  //   std::vector <int> nCombi(npp);
  //   for( int i=0; i<npp; ++i ){ 
  //     nCombi[i]=(CandCont[i]).size();
  //     // If #Cluster>MaxNumberOfCluster,  error return
  //     if( nCombi[i]>MaxNumberOfClusters ){
  //       for( int i=0; i<npp; ++i )
  //   	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  //       return 0;
  //     } 
  //   }
  
  //  if( nCombi[2]==0 ){
#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<npp; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<npp; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
      //      std::cout << CandCont[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif
  //       }
  std::vector < std::vector <int> > 
    CombiIndex = makeindex( npp, &nCombi[0] );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi ) return 0;

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]));
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){
      TrackCont.push_back(track);
    }    
    else
      delete track;
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks

  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<npp; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

//For SDC3&4
// int SdcOutLocalTrackSearch( const DCHitContainer * HC,
// 			    std::vector <DCLocalTrack *> &TrackCont)
		      
// {
//   static const std::string funcname = "[LocalTrackSearch]";

//   std::vector <int> UsedPlane;
//   std::vector <int> UnUsedPlane;

//   double meanMultiSdc3=0.0;
//   for (int i=1; i<=6; i++) {
//     meanMultiSdc3 += (double)HC[i].size();
//   }
//   meanMultiSdc3 /= 6.;
//   for (int i=1; i<=6; i++) {
//     if ((HC[i].size() < meanMultiSdc3+2.) && UsedPlane.size()<4)
//       UsedPlane.push_back(i);
//     else
//       UnUsedPlane.push_back(i);
//   }

//   double meanMultiSdc4=0.0;
//   for (int i=7; i<=12; i++) {
//     meanMultiSdc4 += (double)HC[i].size();
//   }
//   meanMultiSdc4 /= 6.;
//   for (int i=12; i>=7; i--) {
//     if ((HC[i].size() < meanMultiSdc4+2.) && (HC[i].size() > 0) && UsedPlane.size()<8)
//       UsedPlane.push_back(i);
//     else
//       UnUsedPlane.push_back(i);
//   }
  
//   std::vector < std::vector <DCPairHitCluster *> > CandCont;
//   CandCont.resize(UsedPlane.size());

//   for( int i=0; i<UsedPlane.size(); ++i ){
//     MakeUnPairPlaneHitCluster( HC[UsedPlane[i]], CandCont[i] );
//   }

//   std::vector <int> nCombi(UsedPlane.size());
//   for( int i=0; i<UsedPlane.size(); ++i ){ 
//     nCombi[i]=(CandCont[i]).size();

//     // If #Cluster>MaxNumerOfCluster,  error return

//     if(nCombi[i]>MaxNumberOfClusters){
//       for( int i=0; i<UsedPlane.size(); ++i )
// 	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
//       return 0;
//     } 
//   }

// #if 0
//   std::cout << funcname << ": Used Plane" << std::endl;
//   for( int i=0; i<UsedPlane.size(); ++i ) std::cout << std::setw(4) << UsedPlane[i];
//   std::cout << std::endl;

//   std::cout << funcname << ": #Hits of each group" << std::endl;
//   for( int i=0; i<UsedPlane.size(); ++i ) std::cout << std::setw(4) << nCombi[i];
//   std::cout << std::endl;
//   for( int i=0; i<UsedPlane.size(); ++i ){
//     int n=CandCont[i].size();
//     std::cout << "[" << std::setw(3) << i << "]: "
// 	      << std::setw(3) << n << " ";
//     for( int j=0; j<n; ++j ){
//       std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
//     }
//     std::cout << std::endl;
//   }
// #endif

//   std::vector < std::vector <int> > 
//     CombiIndex = makeindex_SdcOut( UsedPlane.size(), 
// 				   MinNumOfHitsSdcOut, 
// 				   UsedPlane.size(), 
// 				   &(nCombi[0]) );
//   int nnCombi=CombiIndex.size();

// #if 0
//   std::cout << " ===> " << nnCombi << " combinations will be checked.." 
// 	    << std::endl;
// #endif

//   if( nnCombi>MaxCombi ) return 0;

//   for( int i=0; i<nnCombi; ++i ){
//     DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
//     if( !track ) continue;
//     if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
// 	track->GetChiSquare()<MaxChisquare ){
//       TrackCont.push_back(track);
//       double chisqr = track->GetChiSquare();
//     }
//     else{
//       //      std::cout << "No tracks available" << std::endl;
//       delete track;
//     }
//   }

//   // Clear Flags
//   int nbefore=TrackCont.size();
//   for( int i=0; i<nbefore; ++i ){
//     DCLocalTrack *tp=TrackCont[i];
//     int nh=tp->GetNHit();
//     for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
//   }

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": Before Sorting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   partial_sort( TrackCont.begin(), TrackCont.end(), 
// 		TrackCont.end(), DCLTrackComp() );

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": After Sorting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   // Delete Duplicated Tracks

//   for( int i=0; i<int(TrackCont.size()); ++i ){
//     DCLocalTrack *tp=TrackCont[i];
//     int nh=tp->GetNHit();
//     for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

//     for( int i2=TrackCont.size()-1; i2>i; --i2 ){
//       DCLocalTrack *tp2=TrackCont[i2];
//       int nh2=tp2->GetNHit(), flag=0;
//       for( int j=0; j<nh2; ++j )
// 	if( tp2->GetHit(j)->showFlags() ) ++flag;
//       if(flag){
// 	delete tp2;
// 	TrackCont.erase(TrackCont.begin()+i2);
//       }
//     }      
//   }

//   {
//     int nn=TrackCont.size();
//     for(int i=0; i<nn; ++i ){
//       DCLocalTrack *tp=TrackCont[i];
//       int nh=tp->GetNHit();
//       for( int j=0; j<nh; ++j ){
// 	int lnum = tp->GetHit(j)->GetLayer();
// 	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
// 	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
//       }
//     }
//   }

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": After Deleting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   int nt=TrackCont.size();
//   for ( int it=0; it<nt; it++ ) {
//     const DCGeomMan & geomMan = DCGeomMan::GetInstance();
//     DCLocalTrack *tp=TrackCont[it];
//     double x0=tp->GetX0(), y0=tp->GetY0();
//     double u0=tp->GetU0(), v0=tp->GetV0();
 
//     for (int ip = 0; ip < UnUsedPlane.size(); ip++) {
//       int lnum = UnUsedPlane[ip]+PlOffsSdcOut;
//       double angle=geomMan.GetTiltAngle( lnum );
//       double z = geomMan.GetLocalZ( lnum );      
//       double x = x0+u0*z;
//       double y = y0+v0*z;
//       double localCalPos = x*cos(angle*Deg2Rad)+y*sin(angle*Deg2Rad);
//       /*
//       std::cout << "Layer# = " << lnum 
// 		<< ", LocalCalPos = " << localCalPos
// 		<< std::endl;
//       */
//       double mindiff = 1000.;
//       int ihUse = -1;
//       double wpUse, dlUse;
//       int mUse;
//       double maxDiff = 5.0;
//       for (int ih=0; ih<HC[UnUsedPlane[ip]].size(); ih++ ){
// 	DCHit *hit = HC[UnUsedPlane[ip]][ih];
// 	if ( hit ) {
// 	  int multi = hit->GetDriftLengthSize();
// 	  for (int m=0; m<multi; m++) {
// 	    if( !(hit->rangecheck(m)) ) continue;
// 	    if( (hit->showFlags(m)) ) continue;
// 	    double wp=hit->GetWirePosition();
// 	    double dl=hit->GetDriftLength(m);
// 	    double diff1 = fabs((wp+dl) - localCalPos);
// 	    double diff2 = fabs((wp-dl) - localCalPos);
// 	    /*
// 	    std::cout << ih << "th Hit " << std::endl;
// 	    std::cout << "localPos1 = " << wp+dl 
// 		      << ", Diff = " << (wp+dl) - localCalPos << std::endl;
// 	    std::cout << "localPos2 = " << wp+dl 
// 		      << ", Diff = " << (wp-dl) - localCalPos << std::endl;
// 	    */
// 	    if (diff1 < maxDiff && diff1 < mindiff ) {
// 	      mindiff = diff1;
// 	      ihUse = ih; mUse = m;
// 	      wpUse = wp; dlUse = dl;
// 	    }
// 	    if (diff2 < maxDiff && diff2 < mindiff ) {
// 	      mindiff = diff2;
// 	      ihUse = ih; mUse = m;
// 	      wpUse = wp; dlUse = -dl;
// 	    }

// 	  }
// 	}
//       }
//       if ( ihUse>=0 && ihUse < HC[UnUsedPlane[ip]].size() ){
// 	DCHit *hit = HC[UnUsedPlane[ip]][ihUse];
// 	DCLTrackHit *lhit = new DCLTrackHit(hit, wpUse+dlUse, mUse);
// 	lhit->setFlags();
// 	if(lhit) tp->AddHit( lhit );
//       }
//     }
//     tp->DoFit();
//   }

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": After Adding UnUsed Plane. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   for( int i=0; i<UsedPlane.size(); ++i )
//     for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
//   return TrackCont.size();
// }

//For SDC3&4
// int SdcOutLocalTrackSearch( const DCHitContainer * HC,
// 			    std::vector <DCLocalTrack *> &TrackCont)
		      
// {
//   static const std::string funcname = "[LocalTrackSearch]";

//   std::vector < std::vector <DCPairHitCluster *> > CandCont;
//   CandCont.resize(NumOfLayersSdcOut);


//   for( int i=0; i<NumOfLayersSdcOut; ++i ){
//     MakeUnPairPlaneHitCluster( HC[i], CandCont[i] );
//   }
  
//   std::vector <int> nCombi(NumOfLayersSdcOut);
//   for( int i=0; i<NumOfLayersSdcOut; ++i ){ 
//     nCombi[i]=(CandCont[i]).size();

//     // If #Cluster>MaxNumerOfCluster,  error return

//     if(nCombi[i]>MaxNumberOfClusters){
//       for( int i=0; i<NumOfLayersSdcOut; ++i )
// 	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
//       return 0;
//     } 
//   }

// #if 0
//   std::cout << funcname << ": #Hits of each group" << std::endl;
//   for( int i=0; i<NumOfLayersSdcOut; ++i ) std::cout << std::setw(4) << nCombi[i];
//   std::cout << std::endl;
//   for( int i=0; i<NumOfLayersSdcOut; ++i ){
//     int n=CandCont[i].size();
//     std::cout << "[" << std::setw(3) << i << "]: "
// 	      << std::setw(3) << n << " ";
//     for( int j=0; j<n; ++j ){
//       std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
//     }
//     std::cout << std::endl;
//   }
// #endif

//   std::vector < std::vector <int> > 
//     CombiIndex = makeindex_SdcOut( NumOfLayersSdcOut, 
// 				   MinNumOfHitsSdcOut, 
// 				   NumOfLayersSdcOut, 
// 				   &(nCombi[0]) );
//   int nnCombi=CombiIndex.size();
 
// #if 0
//   std::cout << " ===> " << nnCombi << " combinations will be checked.." 
// 	    << std::endl;
// #endif

//   for( int i=0; i<nnCombi; ++i ){
//     DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
//     if( !track ) continue;
//     if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
// 	track->GetChiSquare()<MaxChisquare ){
//       TrackCont.push_back(track);
//       double chisqr = track->GetChiSquare();
//     }
//     else{
//       //      std::cout << "No tracks available" << std::endl;
//       delete track;
//     }
//   }

//   // Clear Flags
//   int nbefore=TrackCont.size();
//   for( int i=0; i<nbefore; ++i ){
//     DCLocalTrack *tp=TrackCont[i];
//     int nh=tp->GetNHit();
//     for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
//   }

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": Before Sorting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   partial_sort( TrackCont.begin(), TrackCont.end(), 
// 		TrackCont.end(), DCLTrackComp() );

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": After Sorting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   // Delete Duplicated Tracks

//   for( int i=0; i<int(TrackCont.size()); ++i ){
//     DCLocalTrack *tp=TrackCont[i];
//     int nh=tp->GetNHit();
//     for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

//     for( int i2=TrackCont.size()-1; i2>i; --i2 ){
//       DCLocalTrack *tp2=TrackCont[i2];
//       int nh2=tp2->GetNHit(), flag=0;
//       for( int j=0; j<nh2; ++j )
// 	if( tp2->GetHit(j)->showFlags() ) ++flag;
//       if(flag){
// 	delete tp2;
// 	TrackCont.erase(TrackCont.begin()+i2);
//       }
//     }      
//   }

//   {
//     int nn=TrackCont.size();
//     for(int i=0; i<nn; ++i ){
//       DCLocalTrack *tp=TrackCont[i];
//       int nh=tp->GetNHit();
//       for( int j=0; j<nh; ++j ){
// 	int lnum = tp->GetHit(j)->GetLayer();
// 	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
// 	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
//       }
//     }
//   }

// #if 0
//   {
//     int nn=TrackCont.size();
//     std::cout << funcname << ": After Deleting. #Tracks = " 
// 	      << nn << std::endl;

//     for( int i=0; i<nn; ++i ){
//       DCLocalTrack *track=TrackCont[i];
//       std::cout << std::setw(3) << i << " #Hits="
// 		<< std::setw(2) << track->GetNHit() 
// 		<< " ChiSqr=" << track->GetChiSquare()
// 		<< std::endl;
//     }
//     std::cout << std::endl;

//   }
// #endif

//   for( int i=0; i<NumOfLayersSdcOut; ++i )
//     for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
//   return TrackCont.size();
// }

//For MWPC
int MWPCLocalTrackSearch( const DCHitContainer * HC,
			  std::vector <DCLocalTrack *> &TrackCont)
		      
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayersBcIn);

  for( int i=0; i<NumOfLayersBcIn; ++i ){
    MakeMWPCPairPlaneHitCluster( HC[i], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayersBcIn);
  for( int i=0; i<NumOfLayersBcIn; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayersBcIn; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersBcIn; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersBcIn; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndex = makeindex_BcIn( NumOfLayersBcIn, 
				 MinNumOfHitsBcIn, 
				 NumOfLayersBcIn, 
				 &(nCombi[0]) );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsBcIn && track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){
      TrackCont.push_back(track);
      double chisqr = track->GetChiSquare();
    }
    else{
      //      std::cout << "No tracks available" << std::endl;
      delete track;
    }
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks

  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayersBcIn; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

int LocalTrackSearchBcOutSdcIn( const DCHitContainer * BcHC,  
                                 const DCPairPlaneInfo * BcPpInfo,
                                 const DCHitContainer * SdcHC,  
                                 const DCPairPlaneInfo * SdcPpInfo,
                                 int BcNpp, int SdcNpp,
                                 std::vector <DCLocalTrack *> &TrackCont,
                                 int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearchBcSdc]";

  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(BcNpp+SdcNpp);


  for( int i=0; i<BcNpp; ++i ){
    bool ppFlag=BcPpInfo[i].flag;
    int layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
    if(ppFlag) 
      MakePairPlaneHitCluster( BcHC[layer1], BcHC[layer2], 
                               BcPpInfo[i].CellSize, CandCont[i] );
    else
      MakeUnPairPlaneHitCluster( BcHC[layer1], CandCont[i] );
  }

  for( int i=0; i<SdcNpp; ++i ){
    bool ppFlag=SdcPpInfo[i].flag;
    int layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
    if(ppFlag) 
      MakePairPlaneHitCluster( SdcHC[layer1], SdcHC[layer2], 
                               SdcPpInfo[i].CellSize, CandCont[i+BcNpp] );
    else
      MakeUnPairPlaneHitCluster( SdcHC[layer1], CandCont[i+BcNpp] );
  }

  std::vector <int> nCombi(BcNpp+SdcNpp);

  for( int i=0; i<BcNpp+SdcNpp; ++i ){ 
    nCombi[i]=(CandCont[i]).size();
    int clNum;
    clNum=1;
    if(nCombi[i]!=clNum){
      for( int i=0; i<BcNpp+SdcNpp; ++i )
        for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<BcNpp+SdcNpp; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<BcNpp+SdcNpp; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
              << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << CandCont[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > CombiIndex;
  for (int i=0;i<1;i++) {
    std::vector <int> index;
    index.reserve(BcNpp+SdcNpp);
    for (int j=0;j<BcNpp+SdcNpp;j++) {
      index.push_back(0);
    }
    CombiIndex.push_back(index);
  }
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
            << std::endl;
  for( int i=0; i<nnCombi; ++i ){
    for (int j=0;j<BcNpp+SdcNpp;j++) {
      std::cout << CombiIndex[i][j] << " ";
    }
    std::cout << std::endl;
  }
#endif

  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
        track->GetChiSquare()<MaxChisquare ) 
      TrackCont.push_back(track);
    else
      delete track;
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }


#if 0
 {
   int nn=TrackCont.size();
   std::cout << funcname << ": Before Sorting. #Tracks = " 
	     << nn << std::endl;

   for( int i=0; i<nn; ++i ){
     DCLocalTrack *track=TrackCont[i];
     std::cout << std::setw(3) << i << " #Hits="
	       << std::setw(2) << track->GetNHit() 
	       << " ChiSqr=" << track->GetChiSquare()
	       << std::endl;
   }
   std::cout << std::endl;

 }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
                TrackCont.end(), DCLTrackComp() );

#if 0
 {
   int nn=TrackCont.size();
   std::cout << funcname << ": After Sorting. #Tracks = " 
	     << nn << std::endl;
   /*
    for( int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
                << std::setw(2) << track->GetNHit() 
                << " ChiSqr=" << track->GetChiSquare()
                << std::endl;
    }
    std::cout << std::endl;
   */
 }
#endif

 // Delete Duplicated Tracks

 for( int i=0; i<int(TrackCont.size()); ++i ){
   DCLocalTrack *tp=TrackCont[i];
   int nh=tp->GetNHit();
   for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

   for( int i2=TrackCont.size()-1; i2>i; --i2 ){
     DCLocalTrack *tp2=TrackCont[i2];
     int nh2=tp2->GetNHit(), flag=0;
     for( int j=0; j<nh2; ++j )
       if( tp2->GetHit(j)->showFlags() ) ++flag;
     if(flag){
       delete tp2;
       TrackCont.erase(TrackCont.begin()+i2);
     }
   }      
 }

 {
   double zK18tgt=DCGeomMan::GetInstance().GetLocalZ( 130 );
   int nn=TrackCont.size();
   for(int i=0; i<nn; ++i ){
     DCLocalTrack *tp=TrackCont[i];
     int nh=tp->GetNHit();
     for( int j=0; j<nh; ++j ){
       int lnum = tp->GetHit(j)->GetLayer();
       double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
       if (lnum>=113&&lnum<=124)
	 zz -= zK18tgt;
       tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
     }
   }
 }

#if 0
 {
   int nn=TrackCont.size();
   std::cout << funcname << ": After Deleting. #Tracks = " 
	     << nn << std::endl;

   for( int i=0; i<nn; ++i ){
     DCLocalTrack *track=TrackCont[i];
     std::cout << std::setw(3) << i << " #Hits="
	       << std::setw(2) << track->GetNHit() 
	       << " ChiSqr=" << track->GetChiSquare()
	       << std::endl;
   }
   std::cout << std::endl;

 }
#endif


 for( int i=0; i<SdcNpp+BcNpp; ++i )
   for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
 return TrackCont.size();
}


bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont )
{
  int nh1=HC1.size(), nh2=HC2.size();
  std::vector <int> UsedFlag(nh2,0);   
  
  for( int i1=0; i1<nh1; ++i1 ){
    DCHit *hit1=HC1[i1];

    double wp1=hit1->GetWirePosition();
    bool flag=false;
    for( int i2=0; i2<nh2; ++i2 ){
      DCHit *hit2=HC2[i2];
      double wp2=hit2->GetWirePosition();
      if( fabs(wp1-wp2)<CellSize ){

	int multi1 = hit1->GetDriftLengthSize();
	int multi2 = hit2->GetDriftLengthSize();
	for(int m1=0; m1<multi1; m1++) {
	  if( !(hit1->rangecheck(m1)) ) continue;
	  for(int m2=0; m2<multi2; m2++) {
	    if( !(hit2->rangecheck(m2)) ) continue;
	    double x1,x2;
	    if( wp1<wp2 ){
	      x1=wp1+hit1->GetDriftLength(m1);
	      x2=wp2-hit2->GetDriftLength(m2);
	    }
	    else{
	      x1=wp1-hit1->GetDriftLength(m1);
	      x2=wp2+hit2->GetDriftLength(m2);
	    }
	    
	    Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,x1,m1),
						  new DCLTrackHit(hit2,x2,m2) ) );
	    flag=true; ++UsedFlag[i2];
	  }
	}
      }
    }      
#if 1
    if(!flag){
      int multi1 = hit1->GetDriftLengthSize();
      for (int m1=0; m1<multi1; m1++) {
	if( !(hit1->rangecheck(m1)) ) continue;
	double dl=hit1->GetDriftLength(m1);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1+dl,m1) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit1,wp1-dl,m1) ) );
      }
    }
#endif
  }
#if 1
  for( int i2=0; i2<nh2; ++i2 ){
    if( UsedFlag[i2]==0 ) {
      DCHit *hit2=HC2[i2];
      int multi2 = hit2->GetDriftLengthSize();
      for (int m2=0; m2<multi2; m2++) {
	if( !(hit2->rangecheck(m2)) ) continue;
	
	double wp=hit2->GetWirePosition();
	double dl=hit2->GetDriftLength(m2);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp+dl,m2) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit2,wp-dl,m2) ) );
      }
    }
  }
#endif

  return true;
}

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont )
{
  int nh=HC.size();

  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftLengthSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;

	double wp=hit->GetWirePosition();
	double dl=hit->GetDriftLength(m);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp+dl,m) ) );
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp-dl,m) ) );
      }
    }
  }

  return true;
}

bool MakeMWPCPairPlaneHitCluster( const DCHitContainer & HC,
				  std::vector <DCPairHitCluster *> & Cont )
{
  int nh=HC.size();
  
  for( int i=0; i<nh; ++i ){
    DCHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetDriftTimeSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;
	double wp=hit->GetWirePosition();
	double dl=hit->GetDriftLength(m);
	Cont.push_back( new DCPairHitCluster( new DCLTrackHit(hit,wp,m) ) );
      }
    }
  }

  return true;
}

std::vector< std::vector<int> > makeindex( int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex( ndim-1, index1+1 );

  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      for( int k=0; k<n3; ++k )
	elem.push_back(index2[j][k]);
      index.push_back(elem);
      int size1=index.size();
    }
  }

  return index;
}  

std::vector< std::vector<int> > 
makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_SdcOut( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}  

std::vector< std::vector<int> > 
makeindex_BcIn( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_BcIn( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}  

DCLocalTrack 
* MakeTrack(  const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
	      const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  DCLocalTrack *tp=new DCLocalTrack();
  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    
  int n=CandCont.size();

  for( int i=0; i<n; ++i ){
    int m=combination[i];
    DCPairHitCluster *cluster=0;
    if(m>=0) cluster=CandCont[i][m];
#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif

    if(cluster){
     int m=cluster->NumberOfHits();
      for(int j=0; j<m; ++j ){
	DCLTrackHit *hitp=cluster->GetHit(j);
	if(hitp) tp->AddHit( hitp );
      }
    }
  }

  return tp;
}

//For SDC3&4
//Yonemoto Program
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont)
		      
{
  static const std::string funcname = "[LocalTrackSearch]";
  
  std::vector < std::vector <DCPairHitCluster *> > CandCont;
  CandCont.resize(NumOfLayersSdcOut);



  //Reject X plane
  int Layer_A=-2 ,Layer_B=-2 , Layer_C=-2 ,Layer_D=-2 ,Layer_E=-2 , Layer_F=-2 ;

#if 1
  int W=0 ,NHIT[NumOfLayersSdcOut] ;

  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    int nh ;
    nh = HC[i].size();
    if(nh==0 && W==5){
      Layer_F = i ;
      W++ ;
    }
    if(nh==0 && W==4){
      Layer_E = i ;
      W++ ;
    }

    if(nh==0 && W==3){
      Layer_D = i ;
      W++ ;
    }

    if(nh==0 && W==2){
      Layer_C = i ;
      W++ ;
    }

    if(nh==0 && W==1){
      Layer_B = i ;
      W++ ;
    }

    if(nh==0 && W==0){
      Layer_A = i ;
      W++ ;
    }

    NHIT[i] = nh ;

#if 0
    if(nh==0) std::cout << "#" << i+1 <<"  nh =  " << nh << " W = " << W << std::endl ;
#endif

  }

#if 0
  if(W>0 && W<4){
    for( int i=1;i<NumOfLayersSdcOut+1; ++i ){
      if((i!=Layer_A) && (i!=Layer_A +1) && (i!=Layer_A -1)){
	if((i!=Layer_B) && (i!=Layer_B +1) && (i!=Layer_B -1)){
	  if((i!=Layer_C) && (i!=Layer_C +1) && (i!=Layer_C -1)){
	    if((i!=Layer_D) && (i!=Layer_D +1) && (i!=Layer_D -1)){
	      if(Layer_A!=-2 && Layer_B!=-2 && Layer_C!=-2 && Layer_D!=-2)
		break ; 
	      if(Layer_A!=-2 && Layer_B!=-2 && Layer_C!=-2 && Layer_D==-2)
		Layer_D = i-1 ;
              if(Layer_A!=-2 && Layer_B!=-2 && Layer_C==-2)
		Layer_C = i-1 ;
              if(Layer_A!=-2 && Layer_B==-2)
		Layer_B = i-1 ;
	      if(Layer_A==-2)
		Layer_A = i-1 ;
	    }
	  }
	}
      }
    }
  }
#endif


  if(W==0){
    //Layer_F = 10 ;//SDC4 X2
    //Layer_E = 8  ;//SDC4 U2
    Layer_D = 9  ;//SDC4 V2
    Layer_C = 7  ;//SDC4 X1
    Layer_B = 4  ;//SDC3 X2
    Layer_A = 2  ;//SDC3 U1
  }

#if 1
  if(W>0 && W<4){

    /* Calcurate the average nh (not include nh=0 plane)   */

    int Sum_nh=0 , Num_nh=0 ;
    
    for(int i=0;i<NumOfLayersSdcOut; i++){
      
      int nh ;
      nh = HC[i].size();
      
      if(nh!=0){
	Sum_nh = nh ;
	Num_nh = Num_nh +1 ;
      }
    }

    double Ave_nh = Sum_nh/Num_nh ;



    for(int i=0;i<NumOfLayersSdcOut;i++){
      if((i!=Layer_A -1) && (i!=Layer_A) && (i!=Layer_A +1) && 
	 (i!=Layer_B -1) && (i!=Layer_B) && (i!=Layer_B +1) && 
	 (i!=Layer_C -1) && (i!=Layer_C) && (i!=Layer_C +1) && 
	 (i!=Layer_D -1) && (i!=Layer_D) && (i!=Layer_D +1) && 
	 (i!=Layer_E -1) && (i!=Layer_E) && (i!=Layer_E +1)){
      	int nh ;
	nh = HC[i].size();
	
	/*
	  if(nh>Ave_nh+2 && W==5){
	  Layer_F = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (fourth) = " << Layer_D << std::endl ;
	}
	if(nh>Ave_nh+2 && W==4){
	  Layer_E = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (fourth) = " << Layer_D << std::endl ;
	} 
	*/ 
	if(nh>Ave_nh+2 && W==3){
	  Layer_D = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (thrid) = " << Layer_C << std::endl ;
	}  
	if(nh>Ave_nh+2 && W==2){
	  Layer_C = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (second) = " << Layer_B << std::endl ;
	}  
	if(nh>Ave_nh+2 && W==1){
	  Layer_B = i ;
	  W++ ;
	  //std::cout << "nh>Ave_nh+2 (first) = " << Layer_A << std::endl ;
	}  
      }
    }

    for(int i=0;i<NumOfLayersSdcOut;i++){
      if((i!=Layer_A -1) && (i!=Layer_A) && (i!=Layer_A +1) && 
	 (i!=Layer_B -1) && (i!=Layer_B) && (i!=Layer_B +1) && 
	 (i!=Layer_C -1) && (i!=Layer_C) && (i!=Layer_C +1) && 
	 (i!=Layer_D -1) && (i!=Layer_D) && (i!=Layer_D +1) && 
	 (i!=Layer_E -1) && (i!=Layer_E) && (i!=Layer_E +1)){

	/*
	if(W==5){
	  Layer_F = i ;
	  W++ ;
	}	
	if(W==4){
	  Layer_E = i ;
	  W++ ;
	} 
	*/ 
	if(W==3){
	  Layer_D = i ;
	  W++ ;
	} 
	if(W==2){
	  Layer_C = i ;
	  W++ ;
	}  
	if(W==1){
	  Layer_B = i ;
	  W++ ;
	}  
      }
    }
  }
#endif
  
#if 0
  if(W>=0){
  std::cout << "Layer_A = " << Layer_A+1 << " Layer_B = " << Layer_B+1 
	    << " Layer_C = " << Layer_C+1 << " Layer_D = " << Layer_D+1 
	    << " Layer_E = " << Layer_E+1 << " Layer_F = " << Layer_F+1 << std::endl;
  }
#endif 

#endif
 
  for( int i=0; i<NumOfLayersSdcOut; ++i ){//NumOfLayersSdcOut=12

    if((i!=Layer_A) && (i!=Layer_B) && (i!=Layer_C) && (i!=Layer_D) && (i!=Layer_E) && (i!=Layer_F))
      MakeUnPairPlaneHitCluster( HC[i], CandCont[i] );//Go to bottom
    
  }
  
  std::vector <int> nCombi(NumOfLayersSdcOut);
  for( int i=0; i<NumOfLayersSdcOut; ++i ){ 
    nCombi[i]=(CandCont[i]).size();//always nCombi[i]=2 (left & right by 1 hit)
    
#if 0
    std::cout << "i=  " << i << "   nCombi = " << nCombi[i] << std::endl;
#endif 
    
    if(nCombi[i]>MaxNumberOfClusters){//MaxNumberOfClusters=10
      for( int i=0; i<NumOfLayersSdcOut; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }
  
  
  /*--------HitWire & Hit Local position of each plane----*/
#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayersSdcOut; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " "
                << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetLocalHitPos() << "  :  " ;
    }
    std::cout << std::endl;
  }
#endif
  
  std::vector < std::vector <int> > 
    CombiIndex = makeindex_SdcOut( NumOfLayersSdcOut, 
				   MinNumOfHitsSdcOut, 
				   NumOfLayersSdcOut, 
				   &(nCombi[0]) );//go to bottom
  int nnCombi=CombiIndex.size();
  
#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  
  for( int i=0; i<nnCombi; ++i ){
    DCLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    
    if( track->GetNHit()>=MinNumOfHitsSdcOut && track->DoFit() &&
    	track->GetChiSquare()<MaxChisquare ){
      
      TrackCont.push_back(track);
      double chisqr = track->GetChiSquare();
      
    }
    else{
      delete track;
    }
  }
  
  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
  }
  
  /*-------------------------------------------------------*/


#if 1 //Add residual plane!!
  {
    int nn=TrackCont.size();    
    
    for( int ii=0; ii<nn; ++ii ){
      DCLocalTrack *track=TrackCont[ii];

#if 0
      std::cout << std::setw(3) << "nn = " << nn << " ii = " << ii << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
#endif
      
      
      for( int i=0; i<NumOfLayersSdcOut; ++i ){
	
	if((i==Layer_A) || (i==Layer_B) || (i==Layer_C)  || (i==Layer_D) || (i==Layer_E) || (i==Layer_F)){
	  int nh ;
	  nh = HC[i].size();
	  double Del_A=0 ,Del_B=0 ,Del_1[nh],Del_2[nh] ,Wp_1[nh],Dl_1[nh],Wp_2[nh],Dl_2[nh] ;
	  
	  int I1[nh],I2[nh] ,IA=0 ,IB=0 ,II;

	  DCLTrackHit *hitA[nh] ;
	  DCLTrackHit *hitB[nh] ;

	  for(int j=0;j<nh;++j ){
 	    DCHit * hit = HC[i][j];

#if 0
	    std::cout << "nh = " << nh << " ------j(nh) = "<< j << " ------" << std::endl ;
#endif

 	    if((hit) && (nh!=0)){

	      int multi = hit->GetDriftLengthSize();
	      
	      double wp[nh][multi] , dl[nh][multi] ,theta[nh][multi],th[nh][multi] , th1[nh][multi] ;
	      double del_1[nh][multi] , del_2[nh][multi] ;
	      
	      for (int m=0; m<multi; m++) { 
		if( !(hit->rangecheck(m)) ) continue;
		
		wp[j][m]=hit->GetWirePosition();
		dl[j][m]=hit->GetDriftLength(m);
                theta[j][m]=hit->GetTiltAngle();
		
		double Z = DCGeomMan::GetInstance().GetLocalZ(31+i);
		
                th[j][m] = (theta[j][m]/180)*3.14159 ; 
		
#if 0 
		if(dl[j][m]<0.001){
		  std::cout << "i = " << i << " nh = " << nh << " j = " << j << " multi = " << multi 
			    << " m = " << m << std::endl;
		  std::cout << " wp = "  
			    << (wp[j][m]+dl[j][m])               
			    << " dl = " 
			    << dl[j][m]  << std::endl ;
		  std::cout << "(wp+dl) = "  
			    << (wp[j][m]+dl[j][m])               
			    << " (wp-dl) = " 
			    << (wp[j][m]-dl[j][m])  << std::endl ; 
		  std::cout << " |X - (wp+dl)| = "  
			    << fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]+dl[j][m]))               
			    << " |X - (wp-dl)| = " 
			    << fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]-dl[j][m]))
			    << std::endl ;
		}
#endif
		
		if((fabs(wp[j][m])>pow(10,-10))){

		del_1[j][m] = fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]+dl[j][m]));   
		del_2[j][m] = fabs(track->GetX(Z)*cos(th[j][m])+track->GetY(Z)*sin(th[j][m]) - (wp[j][m]-dl[j][m]));
	      }


	      else{

		del_1[j][m] = 0.0;   
		del_2[j][m] = 0.0;
	      }
	      } 
	      
	      
	      for(int m=0; m<multi; m++){
		if( (hit->rangecheck(m)) ) continue;

		wp[j][m] = 0.0 ; 
		dl[j][m] = 0.0 ;
		del_1[j][m] = 0.0 ;
		del_2[j][m] = 0.0 ;

#if 0
		std::cout << "hit->rangecheck(m)" << std::endl ;
#endif

	      }
	      

	      Del_1[j]=del_1[j][0],Del_2[j]=del_2[j][0] ;
	      
              if(multi==1){

		Dl_1[j] = dl[j][0]; 
		Wp_1[j] = wp[j][0]; 
		
		Dl_2[j] = dl[j][0]; 
		Wp_2[j] = wp[j][0]; 
		
		I1[j] = 0 ;
		I2[j] = 0 ;
		
#if 0
		std::cout << "dl[j][0] = " << dl[j][0] << std::endl ;
		std::cout << "wp[j][0] = " << wp[j][0] << std::endl ;
		std::cout << "del_1[j][0] = " << del_1[j][0] << std::endl ;
		std::cout << "del_2[j][0] = " << del_2[j][0] << std::endl ;
		std::cout << "Wp_1[j]+Dl_1[j] = " << Wp_1[j]+Dl_1[j] << std::endl ;
		std::cout << "Wp_2[j]-Dl_2[j] = " << Wp_2[j]-Dl_2[j] << std::endl ;
#endif
		
	      }
	      
	      if(multi>1){

		//Del_1[j]=del_1[j][0],Del_2[j]=del_2[j][0] ;

		for(int m=0; m<multi-1 ;m++){
		  for(int M=m+1 ; M<multi ;M++){
		    if((del_1[j][M]>pow(10,-10)) && (Del_1[j] > del_1[j][M])){
		      Del_1[j] = del_1[j][M] ;
		      I1[j] = M ; 
		    }
		  }
		}		
		Dl_1[j] = dl[j][I1[j]] ;
		Wp_1[j] = wp[j][I1[j]] ;
		
		for(int m=0; m<multi-1 ;m++){
		  for(int M=m+1 ; M<multi ;M++){
		    if((del_2[j][M]>pow(10,-10)) && (Del_2[j] > del_2[j][M])){
		      Del_2[j] = del_2[j][M] ;
		      I2[j] = M ;
		    }
		  }
		}

		Dl_2[j] = dl[j][I2[j]] ;
		Wp_2[j] = wp[j][I2[j]] ;

#if 0
		std::cout << " j = " << j << " Wp_1[j]+Dl_1[j] = " << Wp_1[j]+Dl_1[j] << std::endl ;	
		std::cout << " j = " << j << " Wp_2[j]-Dl_2[j] = " << Wp_2[j]+Dl_2[j] << std::endl ;
#endif

	      }
	    }

	    if(!((hit) && (nh!=0))){

	      Del_1[j] = 0.0 ;
	      Del_2[j] = 0.0 ; 
	      Wp_1[j] = 0.0 ;
	      Dl_1[j] = 0.0 ;
	      Wp_2[j] = 0.0 ;
	      Dl_2[j] = 0.0 ;

#if 0
	      std::cout << "!((hit) && (nh!=0))" << std::endl ;	      
#endif

	    }

	  }

	  if(nh>0){ 
	    for(int I=0; I<nh; I++){
	      if(Del_1[I]>pow(10,-10)){     
		Del_A = Del_1[I];
                IA = I ;
                break ;
	      }
	    }
	    for(int I=0; I<nh; I++){
	      if(Del_2[I]>pow(10,-10)){     
		Del_B = Del_2[I];
                IB = I ;
                break ;
	      }
	    }


#if 0
	    std::cout << "----------------------------" << std::endl ;

	    std::cout << "Del_A = " << Del_A << std::endl ;
	    std::cout << "Del_B = " << Del_B << std::endl ;
#endif


	    if(nh>1){

	      for(int m=0; m<nh-1 ;m++){
		for(int M=m+1 ; M<nh ;M++){
		  if((Del_1[M]>pow(10,-10)) && (Del_A > Del_1[M])){
		    Del_A = Del_1[M] ;
		    IA = M ; 
		  }
		}
	      }
	      
	      for(int m=0; m<nh-1 ;m++){
		for(int M=m+1 ; M<nh ;M++){
		  if((Del_2[M]>pow(10,-10)) && (Del_B > Del_2[M])){
		    Del_B = Del_2[M] ;
		    IB = M ;
		  }
		}
	      }
	    }

#if 0
	    std::cout << "********************************" << std::endl ;
	    std::cout << "IA = " << IA << " IB = " << IB << std::endl ;
#endif


            if((Del_A>pow(10,-10)) && (Del_B>pow(10,-10))){

	      if((Del_A<=Del_B) && (Del_A<2.0)){
		DCHit *Hit_A = HC[i][IA] ;

#if 0
		std::cout << "I1[IA] = " << I1[IA] << "Del_A = " << Del_A << std::endl ;
		std::cout << "Wp_1[IA]+Dl_1[IA] = " << Wp_1[IA]+Dl_1[IA] << std::endl;
		std::cout << "wp = " << Hit_A->GetWirePosition() << " dl = " << Hit_A->GetDriftLength(I1[IA]) << std::endl;
		std::cout << "wp+dl = " << Hit_A->GetWirePosition()+Hit_A->GetDriftLength(I1[IA]) << std::endl; 
#endif

		DCLTrackHit *hit_A = new DCLTrackHit(Hit_A,Wp_1[IA]+Dl_1[IA],I1[IA]);

		track->AddHit(hit_A);
	      }

	      if((Del_A>Del_B) && (Del_B<2.0)){ 
		DCHit *Hit_B = HC[i][IB] ; 

#if 0
		std::cout << "I2[IB] = " << I2[IB] << " Del_B = " << Del_B << std::endl ;
		std::cout << "Wp_2[IB]-Dl_2[IB] = " << Wp_2[IB]-Dl_2[IB] << std::endl ;
		std::cout << "wp = " << Hit_B->GetWirePosition() << " dl = " << Hit_B->GetDriftLength(I2[IB]) << std::endl;
		std::cout << "wp-dl = " << Hit_B->GetWirePosition()-Hit_B->GetDriftLength(I2[IB]) << std::endl;
#endif

		DCLTrackHit *hit_B = new DCLTrackHit(Hit_B,Wp_2[IB]-Dl_2[IB],I2[IB]);
		
		track->AddHit(hit_B);
	      }
	    
	    }
	  }
	}
      }
    }
  
    
    for(int i=0; i<nn; ++i ){
      DCLocalTrack *track=TrackCont[i];
   
      if( !track ) continue;
      track->DoFit();

#if 0
      if((track->GetChiSquare()>0)) std::cout << "i(nn) = " << i 
						<< " Chisquare = " << track->GetChiSquare() 
						<< " N Hit = " << track->GetNHit() << std::endl ;
#endif

    }
    
    {
      for(int i=0; i<nn; ++i ){
	DCLocalTrack *tp=TrackCont[i];
	int nh=tp->GetNHit();

#if 0
	if(nh<12) std::cout << "---------- N Hit = " << nh << " ----------"<< std::endl ;
#endif


	for( int j=0; j<nh; ++j ){
	  int lnum = tp->GetHit(j)->GetLayer();
	  double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	  tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
	}
      }
    }    
  }    
#endif


  /*-------------------------------------------------------*/



#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;
    if(nn>0){

      if(nn>20) nn = 20 ;
      for( int i=0; i<nn; ++i ){
	//for( int i=0; i<20; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), DCLTrackComp() );
  
  
#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;
    if(nn>0){
      if(nn>20) nn = 20 ;
      for( int i=0; i<nn; ++i ){
	DCLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif
  
  // Delete Duplicated Tracks
  
  int nTC = TrackCont.size(); 
  
#if 0
  std::cout << "nTC =  " << nTC << std::endl ;
#endif
  
  
  for( int i=0; i<int(TrackCont.size()); ++i ){
    DCLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      DCLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags()) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }
  
  {
    int nn=TrackCont.size();

#if 0
  std::cout << "nn =  " << nn << std::endl ;
#endif

    for(int i=0; i<nn; ++i ){
      DCLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();

#if 0
      if((tp->GetChiSquare()>0)) std::cout << "i(nn) = " << i 
					     << " Chisquare(after) = " << tp->GetChiSquare() 
					     << " N Hit = " << nh << std::endl ;
#endif

      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = DCGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }
  
  
  for( int i=0; i<NumOfLayersSdcOut; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );

  return TrackCont.size();
}
