/*
  TrTrackSearch.cc

  2012/5  K.Shirotori
*/

#include "TrTrackSearch.hh"
#include "TrParameters.hh"
#include "TrLTrackHit.hh"
#include "TrHitCluster.hh"
#include "TrLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorID.hh"
#include "TrGeomMan.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#define check_TrTrackSearch 0
#define check_SearchEvent 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double MaxChisquare = 100.;//30
const double MaxChisquareTr = 100.;//30
const double MaxNumberOfClusters = 100.;//10.
const double MaxCombi = 1.0e6;

const double MaxChisquareVXU = 100.;
const double ChisquareCutVXU = 100.;
int NhitGroup=0;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearch( const TrHitContainer * HC,
		      std::vector <TrLocalTrack *> &TrackCont, 
		      int NumOfLayers, int MinNumOfHits,
		      bool FlgMinSegs )
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <TrHitCluster *> > CandCont;
  CandCont.resize(NumOfLayers);

  for( int i=0; i<NumOfLayers; ++i ){
    MakeHitCluster( HC[i], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayers);
  for( int i=0; i<NumOfLayers; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayers; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if check_SearchEvent
// #if check_TrTrackSearch
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  int prod = 1;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    if(n!=0){
      prod *= n;
    }//if(n!=0)
  }//for(i:NumOfLayers)

  std::vector < std::vector <int> > CombiIndex;
  int nnCombi = MaxCombi+1;
  if(prod<=MaxCombi){
    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCombi[0]) );
    nnCombi=CombiIndex.size();
  }else{
    return 0;
  }

#if check_SearchEvent
// #if check_TrTrackSearch
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;

  if(!FlgMinSegs){
    for( int i=0; i<nnCombi; ++i ){
      TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && 
	  track->DoFit() &&
	  track->GetChiSquare()<MaxChisquare ){
	TrackCont.push_back(track);
	double chisqr = track->GetChiSquare();
      }else{
#if check_TrTrackSearch
	std::cout << "No tracks available@" << i+1 << std::endl;
	if(track->GetNHit()<MinNumOfHits){
	  std::cout << "  * Too few hits : " << track->GetNHit() << " < " << MinNumOfHits << std::endl;
	}
	if(!track->DoFit()){
	  std::cout << "  * Wrong tracking : " << track->DoFit() <<  std::endl;
	}
	if(track->GetChiSquare()>=MaxChisquare){
	  std::cout << "  * Large Chesquare : " << track->GetChiSquare() << " >= " << MaxChisquare << std::endl;
	}
#endif
	delete track;
      }
    }//for(i:nnCombi)
  }else{
    for( int i=0; i<nnCombi; ++i ){
      TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
      if( !track ) continue;
      if( track->GetNHit()>=MinNumOfHits && 
	  track->DoFitMS() &&
	  track->GetChiSquare()<MaxChisquare ){
	TrackCont.push_back(track);
	double chisqr = track->GetChiSquare();
      }else{
#if check_TrTrackSearch
	std::cout << "No tracks available@" << i+1 << std::endl;
	if(track->GetNHit()<MinNumOfHits){
	  std::cout << "  * Too few hits : " << track->GetNHit() << " < " << MinNumOfHits << std::endl;
	}
	if(!track->DoFitMS()){
	  std::cout << "  * Wrong tracking : " << track->DoFitMS() <<  std::endl;
	}
	if(track->GetChiSquare()>=MaxChisquare){
	  std::cout << "  * Large Chesquare : " << track->GetChiSquare() << " >= " << MaxChisquare << std::endl;
	}
#endif
	delete track;
      }
    }//for(i:nnCombi)
  }//if(!FlgLowestSeg)

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit(); 
    for( int j=0; j<nh; ++j ) {
      //if(tp->GetHit(j)->showFlags()) 
	tp->GetHit(j)->clearFlags();

    }
  }

// #if check_SearchEvent
#if check_TrTrackSearch
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp() );

// #if check_SearchEvent
#if check_TrTrackSearch
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
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
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
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
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

// #if check_SearchEvent
#if check_TrTrackSearch
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearchQ( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont,
		       int NumOfLayers, int MinNumOfHits )
		      
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <TrHitCluster *> > CandCont;
  CandCont.resize(NumOfLayers);

  for( int i=0; i<NumOfLayers; ++i ){
    MakeHitCluster( HC[i], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayers);
  for( int i=0; i<NumOfLayers; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayers; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCombi[0]) );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;

  for( int i=0; i<nnCombi; ++i ){
    TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && 
	track->DoFit2() &&
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
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit(); 
    for( int j=0; j<nh; ++j ) {
      //if(tp->GetHit(j)->showFlags()) 
	tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
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
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
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
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
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
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearchH( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont,
		       int NumOfLayers, int MinNumOfHits )
		      
{
  static const std::string funcname = "[LocalTrackSearch]";

  std::vector < std::vector <TrHitCluster *> > CandCont;
  CandCont.resize(NumOfLayers);

  for( int i=0; i<NumOfLayers; ++i ){
    MakeHitCluster( HC[i], CandCont[i] );
  }
  
  std::vector <int> nCombi(NumOfLayers);
  for( int i=0; i<NumOfLayers; ++i ){ 
    nCombi[i]=(CandCont[i]).size();

    // If #Cluster>MaxNumerOfCluster,  error return

    if(nCombi[i]>MaxNumberOfClusters){
      for( int i=0; i<NumOfLayers; ++i )
	for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > 
    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCombi[0]) );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;

  for( int i=0; i<nnCombi; ++i ){
    TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && 
	track->DoFit4() &&
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
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit(); 
    for( int j=0; j<nh; ++j ) {
      //if(tp->GetHit(j)->showFlags()) 
	tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
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
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
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
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
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
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}

int LocalTrackSearch2( const TrHitContainer * HC,
		      std::vector <TrLocalTrack *> &TrackCont,
		      int NumOfLayers, int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  int MinNumOfHitsVXU = 3;
  std::vector <TrLocalTrack *>  TrackContV;
  std::vector <TrLocalTrack *>  TrackContX;
  std::vector <TrLocalTrack *>  TrackContU;

  std::vector < std::vector <TrHitCluster *> > CandCont;
  std::vector < std::vector <TrHitCluster *> > CandContV;
  std::vector < std::vector <TrHitCluster *> > CandContX;
  std::vector < std::vector <TrHitCluster *> > CandContU;

  CandCont.resize(NumOfLayers);
  CandContV.resize(NumOfLayers/3);
  CandContX.resize(NumOfLayers/3);
  CandContU.resize(NumOfLayers/3);

  for( int i=0; i<NumOfLayers/3; ++i ){
    MakeHitCluster( HC[3*i+1], CandContV[i] );
    MakeHitCluster( HC[3*i+2], CandContX[i] );
    MakeHitCluster( HC[3*i+3], CandContU[i] );
  }   

  std::vector <int> nCombi(NumOfLayers);
  std::vector <int> nCombiV(NumOfLayers/3);
  std::vector <int> nCombiX(NumOfLayers/3);
  std::vector <int> nCombiU(NumOfLayers/3);

  int nV=0, nX=0, nU=0;

  for( int i=0; i<NumOfLayers/3; ++i ){ 
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();
    if(nCombiV[i]>0) ++nV;
    if(nCombiX[i]>0) ++nX;
    if(nCombiU[i]>0) ++nU;

    if(nCombiV[i]>MaxNumberOfClusters || nCombiX[i]>MaxNumberOfClusters || nCombiU[i]>MaxNumberOfClusters){
      for( int ii=0; ii<NumOfLayers; ++ii ){
	for_each( CandContV[ii].begin(), CandContV[ii].end(), DeleteObject() );
	for_each( CandContX[ii].begin(), CandContX[ii].end(), DeleteObject() );
	for_each( CandContU[ii].begin(), CandContU[ii].end(), DeleteObject() );
	return 0;
      }
    } 
  }

  //  const  int NhitGroup=0;
  for(int i=0; i<NumOfLayers/3; ++i){
    if(nCombiV[i]>3 || nCombiX[i]>3 || nCombiU[i]>3 ){
      NhitGroup++;
    }
  }

#if 0
  //////////////////----------------  V Plane -------------------////////////////////////////
  std::cout << funcname << ": V plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiV[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContV[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContV[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  X Plane -------------------////////////////////////////
  std::cout << funcname << ": X plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiX[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContX[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContX[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  U Plane -------------------////////////////////////////
  std::cout << funcname << ": U plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiU[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContU[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContU[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > CombiIndexV;
  if(nV>=MinNumOfHitsVXU){
    CombiIndexV = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiV[0]) );
  }

  std::vector < std::vector <int> > CombiIndexX;
  if(nX>=MinNumOfHitsVXU){
    CombiIndexX = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiX[0]) );
  }

  std::vector < std::vector <int> > CombiIndexU;
  if(nU>=MinNumOfHitsVXU){
    CombiIndexU = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiU[0]) );
  }

  int nnCombiV=0, nnCombiX=0, nnCombiU=0;

  if(nV>=MinNumOfHitsVXU) nnCombiV=CombiIndexV.size();
  if(nX>=MinNumOfHitsVXU) nnCombiX=CombiIndexX.size();
  if(nU>=MinNumOfHitsVXU) nnCombiU=CombiIndexU.size();
  
#if 0
  std::cout << "V Plane  ===> " << nnCombiV << " combinations will be checked.." 
	    << std::endl;
  std::cout << "X Plane  ===> " << nnCombiX << " combinations will be checked.." 
	    << std::endl;
  std::cout << "U Plane  ===> " << nnCombiU << " combinations will be checked.." 
	    << std::endl;
#endif

  int NHitV=0,NHitX=0,NHitU=0,NResV=0,NResX=0,NResU=0;

  //for( int i=0; i<1; ++i ){
  for( int i=0; i<nnCombiV; ++i ){
    TrLocalTrack *track = MakeTrack( CandContV, &((CombiIndexV[i])[0]) );

    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContV.push_back(track);
      NHitV = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  for( int i=0; i<nnCombiX; ++i ){
    TrLocalTrack *track = MakeTrack( CandContX, &((CombiIndexX[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContX.push_back(track);
      NHitX = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  for( int i=0; i<nnCombiU; ++i ){
    TrLocalTrack *track = MakeTrack( CandContU, &((CombiIndexU[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContU.push_back(track);
      NHitU = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  
  int nbeforeV=0, nbeforeX=0, nbeforeU=0;
    // Clear Flags
  if(nV>=MinNumOfHitsVXU){
    nbeforeV=TrackContV.size();
    for( int i=0; i<nbeforeV; ++i ){
      TrLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nX>=MinNumOfHitsVXU){
    nbeforeX=TrackContX.size();
    for( int i=0; i<nbeforeX; ++i ){
      TrLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nU>=MinNumOfHitsVXU){
    nbeforeU=TrackContU.size();
    for( int i=0; i<nbeforeU; ++i ){
      TrLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] Before Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] Before Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] Before Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  //  if(nbeforeV>0)  
  partial_sort( TrackContV.begin(), TrackContV.end(), 
		TrackContV.end(), TrLTrackComp1() );
  //  if(nbeforeX>0)  
  partial_sort( TrackContX.begin(), TrackContX.end(), 
		TrackContX.end(), TrLTrackComp1() );
  //  if(nbeforeU>0)  
  partial_sort( TrackContU.begin(), TrackContU.end(), 
		TrackContU.end(), TrLTrackComp1() );

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Av)=" << track->GetVXU_A()*cos(Deg2Rad*(-30.0))
		  << " dY/dZ(Av)=" << (-1)*track->GetVXU_A()*sin(Deg2Rad*(-30.0))
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Ax)=" << track->GetVXU_A()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Au)=" << track->GetVXU_A()*cos(Deg2Rad*30.0)
		  << " dY/dZ(Au)=" << track->GetVXU_A()*sin(Deg2Rad*30.0)
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  double chiV=ChisquareCutVXU, chiX=ChisquareCutVXU, chiU=ChisquareCutVXU;

#if 1
  {
    /* V Plane  */  
    for( int i=0; i<int(TrackContV.size()); ++i ){
      TrLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      if(chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      //      if(nh==nV && chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      else{
	delete tp;
	TrackContV.erase(TrackContV.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContV.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContV[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContV.erase(TrackContV.begin()+i2);
	}
      }      
    }
    /* X Plane  */    
    for( int i=0; i<int(TrackContX.size()); ++i ){
      TrLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      if(chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      //      if(nh==nX && chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      else{
	delete tp;
	TrackContX.erase(TrackContX.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContX.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContX[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContX.erase(TrackContX.begin()+i2);
	}
      }      
    }
    /* U Plane  */    
    for( int i=0; i<int(TrackContU.size()); ++i ){
      TrLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      if(chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      //      if(nh==nU && chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      else{
	delete tp;
	TrackContU.erase(TrackContU.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContU.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContU[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContU.erase(TrackContU.begin()+i2);
	}
      }      
    }
  }
#endif


#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Delete. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Delete. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;

      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Delete. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  int nnV=1,nnX=1,nnU=1;
  int nnVT=1,nnXT=1,nnUT=1;
  int checkV=0,checkX=0,checkU=0;

  int nkV=0,nkX=0,nkU=0;

  int cV=TrackContV.size();
  if(chiV>1.5 || cV<5. ) checkV++;
  int cX=TrackContX.size();
  if(chiX>1.5 || cX<5. ) checkX++;
  int cU=TrackContU.size();
  if(chiU>1.5 || cU<5. ) checkU++;

  std::vector < std::vector <int> > CombiIndexSV;
  std::vector < std::vector <int> > CombiIndexSX;
  std::vector < std::vector <int> > CombiIndexSU;

  {
    if((nV>=MinNumOfHitsVXU) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      //      CombiIndexSV = makeindex_below( NumOfLayers/3, 
      CombiIndexSV = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiV[0]) );
      nnV = nnV +  CombiIndexSV.size();
    }
  }

  {
    if((nX>=MinNumOfHitsVXU) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiX[0]) );
      nnX = nnX + CombiIndexSX.size();
    }
  }

  {
    if((nU>=MinNumOfHitsVXU) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiU[0]) );
      nnU = nnU + CombiIndexSU.size();
    }
  }

  /*** Event Cut by using anular difference***/
  double DifVXU=0.0;
  double Av=0.0;
  double Ax=0.0;
  double Au=0.0;

  bool angularVcut[nnVT+1];
  bool angularXcut[nnXT+1];
  bool angularUcut[nnUT+1];


  for( int i=0; i<nnVT; ++i){
    for( int j=0; j<nnXT; ++j){
      for( int k=0; k<nnUT; ++k){

	angularVcut[i]=false;
	angularXcut[j]=false;
	angularUcut[k]=false;

	if(nkV){
	  TrLocalTrack *trackV=TrackContV[i];
	  Av=trackV->GetVXU_A();
	}
	if(nkX){
	  TrLocalTrack *trackX=TrackContX[j];
	  Ax=trackX->GetVXU_A();
	}
	if(nkU){
	  TrLocalTrack *trackU=TrackContU[k];
	  Au=trackU->GetVXU_A();
	}

	DifVXU = (Av*cos(acos(-1.)/180.*(-30.0))-Ax)*(Av*cos(acos(-1.)/180.*(-30.0))-Ax)
	  +(Ax-Au*cos(acos(-1.)/180.*(30.0)))*(Ax-Au*cos(acos(-1.)/180.*(30.0)))
	  +(Au*cos(acos(-1.)/180.*(30.0))
	    -Av*cos(acos(-1.)/180.*(-30.0)))*(Au*cos(acos(-1.)/180.*(30.0))
					      -Av*cos(acos(-1.)/180.*(-30.0)))
	  + (Au*sin(acos(-1.)/180.*(30.0))
	     -Av*sin(acos(-1.)/180.*(30.0)))*(Au*sin(acos(-1.)/180.*(30.0))
					      -Av*sin(acos(-1.)/180.*(30.0)));

	if( (DifVXU>=0) ){
	  angularVcut[i]=true;
	  angularXcut[j]=true;
	  angularUcut[k]=true;
	}
      }
    }
  }

  DifVXU=0.0;
  Av=0.0;
  Ax=0.0;
  Au=0.0;

  /////////////////////////////
  //V,X,U Tracks -> Add Track// 
  /////////////////////////////

 double chiv, chix, chiu;

  for( int i=-1; i<nnV; ++i){
    for( int j=-1; j<nnX; ++j){
      for( int k=-1; k<nnU; ++k){
	chiv=-1.0,chix=-1.0,chiu=-1.0;
	TrLocalTrack *track = new TrLocalTrack();
	
	int mV=0,mX=0,mU=0;
	int mmV=0,mmX=0,mmU=0;
	
	/* V Plane  */
	if(i>-1){
	  if(nkV && i<nnVT){
	    TrLocalTrack *trackV=TrackContV[i];
	    Av=trackV->GetVXU_A();
	    chiv=trackV->GetChiSquare();
	    
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      TrLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	  }
	  if((i>=nnVT) && (nV>0)){
	    TrLocalTrack *trackV = MakeTrack( CandContV, &((CombiIndexSV[i-nnVT])[0]) ); 
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      TrLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    delete trackV ;
	  }
	  int NHitV = track->GetNHit();
	}

	/* X Plane  */
	if(j>-1){
	  if(nkX && j<nnXT){
	    TrLocalTrack *trackX=TrackContX[j];
	    Ax=trackX->GetVXU_A();
	    chix=trackX->GetChiSquare();
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      TrLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	  }
	  if((j>=nnXT) && (nX>0)){
	    TrLocalTrack *trackX = MakeTrack( CandContX, &((CombiIndexSX[j-nnXT])[0]) ); 
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      TrLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    delete trackX;
	  }


	  int NHitX = track->GetNHit();
	}
	/* U Plane  */
	if(k>-1){
	  if(nkU && k<nnUT){
	    TrLocalTrack *trackU=TrackContU[k];
	    Au=trackU->GetVXU_A();
	    chiu=trackU->GetChiSquare();
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      TrLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	  }
	  if((k>=nnUT) && (nU>0)){
	    TrLocalTrack *trackU = MakeTrack( CandContU, &((CombiIndexSU[k-nnUT])[0]) ); 
	    
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      TrLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    delete trackU;
	  }
	  int NHitU = track->GetNHit();
	}
	
	track->SetAv(Av);
	track->SetAx(Ax);
	track->SetAu(Au);
	DifVXU = track->GetDifVXU2();
	
	track->SetChiv(chiv);
	track->SetChix(chix);
	track->SetChiu(chiu);
	
	if(!track) continue;
	if(track->GetNHit()>=MinNumOfHits && 
	   track->DoFit() && 
	   track->GetChiSquare()<MaxChisquareTr ){//MaXChisquare
	  
	  TrackCont.push_back(track);
	}    
	else{
	  delete track;
	}
      }
    }
  }
  
  for( int i=0; i<NumOfLayers/3; ++i ){
    for_each( CandContV[i].begin(), CandContV[i].end(), DeleteObject() );
    for_each( CandContX[i].begin(), CandContX[i].end(), DeleteObject() );
    for_each( CandContU[i].begin(), CandContU[i].end(), DeleteObject() );
  }

  {
    for_each( TrackContV.begin(), TrackContV.end(), DeleteObject() );
    for_each( TrackContX.begin(), TrackContX.end(), DeleteObject() );
    for_each( TrackContU.begin(), TrackContU.end(), DeleteObject() );
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp1() );    

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " DifVXU=" << track->GetDifVXU2()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

#if 1  
  // Delete Tracks about (Nhit1 = Nhit2  && chi1 < chi2)
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
      double chi=tp->GetChiSquare();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      double chi2=tp2->GetChiSquare();
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if((flag) && ((nh==nh2) && (chi<=chi2))){
        delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }    
#endif  

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": Before Deleted Sorting. #Tracks = " 
		<< nn << std::endl;
      
      //for( int i=0; i<100; ++i ){
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackCont[i];
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
  		TrackCont.end(), TrLTrackComp() );

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
      		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
      //for( int i=0; i<100; ++i ){
	TrLocalTrack *track=TrackCont[i];
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
  
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;  
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }  

  
  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }
  
#if 0
  {
    int nn=TrackCont.size();
    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      if((track->GetChiSquare()>30)){
	// std::cout << "" << std::endl;
	// std::cout << std::setw(3) <<"ntrack = " << nn << " Track# =  " << i << " #Hits="
	// 	  << std::setw(2) << track->GetNHit() 
	// 	  << " ChiSqr=" << track->GetChiSquare()
	// 	  << " DifVXU=" << track->GetDifVXU()
	// 	  << std::endl;
	// std::cout << "" << std::endl;
      }
    }
  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );


  return TrackCont.size();    
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearchQ2( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont,
		       int NumOfLayers, int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";

  int MinNumOfHitsVXU = 7;
  std::vector <TrLocalTrack *>  TrackContV;
  std::vector <TrLocalTrack *>  TrackContX;
  std::vector <TrLocalTrack *>  TrackContU;

  std::vector < std::vector <TrHitCluster *> > CandCont;
  std::vector < std::vector <TrHitCluster *> > CandContV;
  std::vector < std::vector <TrHitCluster *> > CandContX;
  std::vector < std::vector <TrHitCluster *> > CandContU;

  CandCont.resize(NumOfLayers);
  CandContV.resize(NumOfLayers/3);
  CandContX.resize(NumOfLayers/3);
  CandContU.resize(NumOfLayers/3);

  for( int i=0; i<NumOfLayers/3; ++i ){
    MakeHitCluster( HC[3*i+1], CandContV[i] );
    MakeHitCluster( HC[3*i+2], CandContX[i] );
    MakeHitCluster( HC[3*i+3], CandContU[i] );
  }   

  std::vector <int> nCombi(NumOfLayers);
  std::vector <int> nCombiV(NumOfLayers/3);
  std::vector <int> nCombiX(NumOfLayers/3);
  std::vector <int> nCombiU(NumOfLayers/3);

  int nV=0, nX=0, nU=0;

  for( int i=0; i<NumOfLayers/3; ++i ){ 
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();
    if(nCombiV[i]>0) ++nV;
    if(nCombiX[i]>0) ++nX;
    if(nCombiU[i]>0) ++nU;

    if(nCombiV[i]>MaxNumberOfClusters || nCombiX[i]>MaxNumberOfClusters || nCombiU[i]>MaxNumberOfClusters){
      for( int ii=0; ii<NumOfLayers; ++ii ){
	for_each( CandContV[ii].begin(), CandContV[ii].end(), DeleteObject() );
	for_each( CandContX[ii].begin(), CandContX[ii].end(), DeleteObject() );
	for_each( CandContU[ii].begin(), CandContU[ii].end(), DeleteObject() );
	return 0;
      }
    } 
  }

  //  const  int NhitGroup=0;
  for(int i=0; i<NumOfLayers/3; ++i){
    if(nCombiV[i]>3 || nCombiX[i]>3 || nCombiU[i]>3 ){
      NhitGroup++;
    }
  }

#if 0
  //////////////////----------------  V Plane -------------------////////////////////////////
  std::cout << funcname << ": V plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiV[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContV[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContV[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  X Plane -------------------////////////////////////////
  std::cout << funcname << ": X plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiX[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContX[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContX[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
  
  //////////////////----------------  U Plane -------------------////////////////////////////
  std::cout << funcname << ": U plane #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ) std::cout << std::setw(4) << nCombiU[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers/3; ++i ){
    int n=CandContU[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandContU[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif

  std::vector < std::vector <int> > CombiIndexV;
  if(nV>=MinNumOfHitsVXU){
    CombiIndexV = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiV[0]) );
  }

  std::vector < std::vector <int> > CombiIndexX;
  if(nX>=MinNumOfHitsVXU){
    CombiIndexX = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiX[0]) );
  }

  std::vector < std::vector <int> > CombiIndexU;
  if(nU>=MinNumOfHitsVXU){
    CombiIndexU = makeindex( NumOfLayers/3, 
				    MinNumOfHitsVXU, 
				    NumOfLayers/3, 
				    &(nCombiU[0]) );
  }

  int nnCombiV=0, nnCombiX=0, nnCombiU=0;

  if(nV>=MinNumOfHitsVXU) nnCombiV=CombiIndexV.size();
  if(nX>=MinNumOfHitsVXU) nnCombiX=CombiIndexX.size();
  if(nU>=MinNumOfHitsVXU) nnCombiU=CombiIndexU.size();
  
#if 0
  std::cout << "V Plane  ===> " << nnCombiV << " combinations will be checked.." 
	    << std::endl;
  std::cout << "X Plane  ===> " << nnCombiX << " combinations will be checked.." 
	    << std::endl;
  std::cout << "U Plane  ===> " << nnCombiU << " combinations will be checked.." 
	    << std::endl;
#endif

  int NHitV=0,NHitX=0,NHitU=0,NResV=0,NResX=0,NResU=0;

  //for( int i=0; i<1; ++i ){
  for( int i=0; i<nnCombiV; ++i ){
    TrLocalTrack *track = MakeTrack( CandContV, &((CombiIndexV[i])[0]) );

    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU2() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContV.push_back(track);
      NHitV = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  for( int i=0; i<nnCombiX; ++i ){
    TrLocalTrack *track = MakeTrack( CandContX, &((CombiIndexX[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU2() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContX.push_back(track);
      NHitX = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  for( int i=0; i<nnCombiU; ++i ){
    TrLocalTrack *track = MakeTrack( CandContU, &((CombiIndexU[i])[0]) );
    
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHitsVXU && track->DoFitVXU2() &&
	track->GetChiSquare()<MaxChisquareVXU ){
      TrackContU.push_back(track);
      NHitU = track->GetNHit();
      double chisqr = track->GetChiSquare();
    }
    else{
      delete track;
    }
  }
  
  int nbeforeV=0, nbeforeX=0, nbeforeU=0;
    // Clear Flags
  if(nV>=MinNumOfHitsVXU){
    nbeforeV=TrackContV.size();
    for( int i=0; i<nbeforeV; ++i ){
      TrLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nX>=MinNumOfHitsVXU){
    nbeforeX=TrackContX.size();
    for( int i=0; i<nbeforeX; ++i ){
      TrLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

  if(nU>=MinNumOfHitsVXU){
    nbeforeU=TrackContU.size();
    for( int i=0; i<nbeforeU; ++i ){
      TrLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] Before Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] Before Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] Before Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  //  if(nbeforeV>0)  
  partial_sort( TrackContV.begin(), TrackContV.end(), 
		TrackContV.end(), TrLTrackComp1() );
  //  if(nbeforeX>0)  
  partial_sort( TrackContX.begin(), TrackContX.end(), 
		TrackContX.end(), TrLTrackComp1() );
  //  if(nbeforeU>0)  
  partial_sort( TrackContU.begin(), TrackContU.end(), 
		TrackContU.end(), TrLTrackComp1() );

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Sorting. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Av)=" << track->GetVXU_A()*cos(Deg2Rad*(-30.0))
		  << " dY/dZ(Av)=" << (-1)*track->GetVXU_A()*sin(Deg2Rad*(-30.0))
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Sorting. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Ax)=" << track->GetVXU_A()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Sorting. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " dX/dZ(Au)=" << track->GetVXU_A()*cos(Deg2Rad*30.0)
		  << " dY/dZ(Au)=" << track->GetVXU_A()*sin(Deg2Rad*30.0)
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  double chiV=ChisquareCutVXU, chiX=ChisquareCutVXU, chiU=ChisquareCutVXU;

#if 1
  {
    /* V Plane  */  
    for( int i=0; i<int(TrackContV.size()); ++i ){
      TrLocalTrack *tp=TrackContV[i];
      int nh=tp->GetNHit();
      if(chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      //      if(nh==nV && chiV > tp->GetChiSquare()) chiV = tp->GetChiSquare();
      else{
	delete tp;
	TrackContV.erase(TrackContV.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContV.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContV[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContV.erase(TrackContV.begin()+i2);
	}
      }      
    }
    /* X Plane  */    
    for( int i=0; i<int(TrackContX.size()); ++i ){
      TrLocalTrack *tp=TrackContX[i];
      int nh=tp->GetNHit();
      if(chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      //      if(nh==nX && chiX > tp->GetChiSquare()) chiX = tp->GetChiSquare();
      else{
	delete tp;
	TrackContX.erase(TrackContX.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContX.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContX[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContX.erase(TrackContX.begin()+i2);
	}
      }      
    }
    /* U Plane  */    
    for( int i=0; i<int(TrackContU.size()); ++i ){
      TrLocalTrack *tp=TrackContU[i];
      int nh=tp->GetNHit();
      if(chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      //      if(nh==nU && chiU > tp->GetChiSquare()) chiU = tp->GetChiSquare();
      else{
	delete tp;
	TrackContU.erase(TrackContU.begin()+i);
	continue;
      }
      for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
      
      for( int i2=TrackContU.size()-1; i2>i; --i2 ){
	TrLocalTrack *tp2=TrackContU[i2];
	int nh2=tp2->GetNHit(), flag=0;
	for( int j=0; j<nh2; ++j )
	  if( tp2->GetHit(j)->showFlags() ) ++flag;
	if((flag) && tp2->GetChiSquare()>ChisquareCutVXU){
	  delete tp2;
	  TrackContU.erase(TrackContU.begin()+i2);
	}
      }      
    }
  }
#endif

#if 0
  {
    /* V Plane  */
    {
      int nn=TrackContV.size();
      std::cout << funcname << ": [V Plane] After Delete. #Tracks = " 
		<< nn << " nV = " << nV << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContV[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
    /* X Plane  */
    {
      int nn=TrackContX.size();
      std::cout << funcname << ": [X Plane] After Delete. #Tracks = " 
		<< nn << " nX = " << nX << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContX[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;

      }
      std::cout << std::endl;
    }
    /* U Plane  */
    {
      int nn=TrackContU.size();
      std::cout << funcname << ": [U Plane] After Delete. #Tracks = " 
		<< nn << " nU = " << nU << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackContU[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

  if( TrackContV.size()>150 && TrackContX.size()>150 && TrackContU.size()>150 ){
    for( int i=0; i<NumOfLayers/3; ++i ){
      for_each( CandContV[i].begin(), CandContV[i].end(), DeleteObject() );
      for_each( CandContX[i].begin(), CandContX[i].end(), DeleteObject() );
      for_each( CandContU[i].begin(), CandContU[i].end(), DeleteObject() );
    }
    {
      for_each( TrackContV.begin(), TrackContV.end(), DeleteObject() );
      for_each( TrackContX.begin(), TrackContX.end(), DeleteObject() );
      for_each( TrackContU.begin(), TrackContU.end(), DeleteObject() );
    }
    return 0;
  }

  int nnV=1,nnX=1,nnU=1;
  int nnVT=1,nnXT=1,nnUT=1;
  int checkV=0,checkX=0,checkU=0;

  int nkV=0,nkX=0,nkU=0;

  int cV=TrackContV.size();
  if(chiV>1.5 || cV<5. ) checkV++;
  int cX=TrackContX.size();
  if(chiX>1.5 || cX<5. ) checkX++;
  int cU=TrackContU.size();
  if(chiU>1.5 || cU<5. ) checkU++;

  std::vector < std::vector <int> > CombiIndexSV;
  std::vector < std::vector <int> > CombiIndexSX;
  std::vector < std::vector <int> > CombiIndexSU;

  {
    if((nV>=MinNumOfHitsVXU) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      CombiIndexSV = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiV[0]) );
      nnV = nnV +  CombiIndexSV.size();
    }
  }

  {
    if((nX>=MinNumOfHitsVXU) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiX[0]) );
      nnX = nnX + CombiIndexSX.size();
    }
  }

  {
    if((nU>=MinNumOfHitsVXU) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = makeindex_below( NumOfLayers/3, 
				      MinNumOfHitsVXU-1, 
				      NumOfLayers/3, 
				      &(nCombiU[0]) );
      nnU = nnU + CombiIndexSU.size();
    }
  }

  /*** Event Cut by using anular difference***/
  double DifVXU=0.0;
  double Av=0.0;
  double Ax=0.0;
  double Au=0.0;

  bool angularVcut[nnVT+1];
  bool angularXcut[nnXT+1];
  bool angularUcut[nnUT+1];

  for( int i=0; i<nnVT; ++i){
    for( int j=0; j<nnXT; ++j){
      for( int k=0; k<nnUT; ++k){

	angularVcut[i]=false;
	angularXcut[j]=false;
	angularUcut[k]=false;

	if(nkV){
	  TrLocalTrack *trackV=TrackContV[i];
	  Av=trackV->GetVXU_A();
	}
	if(nkX){
	  TrLocalTrack *trackX=TrackContX[j];
	  Ax=trackX->GetVXU_A();
	}
	if(nkU){
	  TrLocalTrack *trackU=TrackContU[k];
	  Au=trackU->GetVXU_A();
	}

	DifVXU = (Av*cos(acos(-1.)/180.*(-30.0))-Ax)*(Av*cos(acos(-1.)/180.*(-30.0))-Ax)
	  +(Ax-Au*cos(acos(-1.)/180.*(30.0)))*(Ax-Au*cos(acos(-1.)/180.*(30.0)))
	  +(Au*cos(acos(-1.)/180.*(30.0))
	    -Av*cos(acos(-1.)/180.*(-30.0)))*(Au*cos(acos(-1.)/180.*(30.0))
					      -Av*cos(acos(-1.)/180.*(-30.0)))
	  + (Au*sin(acos(-1.)/180.*(30.0))
	     -Av*sin(acos(-1.)/180.*(30.0)))*(Au*sin(acos(-1.)/180.*(30.0))
					      -Av*sin(acos(-1.)/180.*(30.0)));

	if( (DifVXU>=0) ){
	  angularVcut[i]=true;
	  angularXcut[j]=true;
	  angularUcut[k]=true;
	}
      }
    }
  }

  DifVXU=0.0;
  Av=0.0;
  Ax=0.0;
  Au=0.0;

  /////////////////////////////
  //V,X,U Tracks -> Add Track// 
  /////////////////////////////

 double chiv, chix, chiu;
  for( int i=-1; i<nnV; ++i){
    for( int j=-1; j<nnX; ++j){
      for( int k=-1; k<nnU; ++k){
	chiv=-1.0,chix=-1.0,chiu=-1.0;
	TrLocalTrack *track = new TrLocalTrack();
	
	int mV=0,mX=0,mU=0;
	int mmV=0,mmX=0,mmU=0;
	
	/* V Plane  */
	if(i>-1){
	  if(nkV && i<nnVT){
	    TrLocalTrack *trackV=TrackContV[i];
	    Av=trackV->GetVXU_A();
	    chiv=trackV->GetChiSquare();
	    
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      TrLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	  }
	  if((i>=nnVT) && (nV>0)){
	    TrLocalTrack *trackV = MakeTrack( CandContV, &((CombiIndexSV[i-nnVT])[0]) ); 
	    for( int l=0; l<(trackV->GetNHit()); ++l){
	      TrLTrackHit *hitpV=trackV->GetHit(l);
	      if(hitpV){
		track->AddHit( hitpV ) ;
	      }
	    }
	    delete trackV ;
	  }
	  int NHitV = track->GetNHit();
	}

	/* X Plane  */
	if(j>-1){
	  if(nkX && j<nnXT){
	    TrLocalTrack *trackX=TrackContX[j];
	    Ax=trackX->GetVXU_A();
	    chix=trackX->GetChiSquare();
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      TrLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	  }
	  if((j>=nnXT) && (nX>0)){
	    TrLocalTrack *trackX = MakeTrack( CandContX, &((CombiIndexSX[j-nnXT])[0]) ); 
	    for( int l=0; l<(trackX->GetNHit()); ++l){
	      TrLTrackHit *hitpX=trackX->GetHit(l);
	      if(hitpX){
		track->AddHit( hitpX ) ;
	      }
	    }
	    delete trackX;
	  }


	  int NHitX = track->GetNHit();
	}
	/* U Plane  */
	if(k>-1){
	  if(nkU && k<nnUT){
	    TrLocalTrack *trackU=TrackContU[k];
	    Au=trackU->GetVXU_A();
	    chiu=trackU->GetChiSquare();
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      TrLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	  }
	  if((k>=nnUT) && (nU>0)){
	    TrLocalTrack *trackU = MakeTrack( CandContU, &((CombiIndexSU[k-nnUT])[0]) ); 
	    
	    for( int l=0; l<(trackU->GetNHit()); ++l){
	      TrLTrackHit *hitpU=trackU->GetHit(l);
	      if(hitpU){
		track->AddHit( hitpU ) ;
	      }
	    }
	    delete trackU;
	  }
	  int NHitU = track->GetNHit();
	}
	
	track->SetAv(Av);
	track->SetAx(Ax);
	track->SetAu(Au);
	DifVXU = track->GetDifVXU2();
	
	track->SetChiv(chiv);
	track->SetChix(chix);
	track->SetChiu(chiu);
	
	if(!track) continue;
	if(track->GetNHit()>=MinNumOfHits && 
	   track->DoFit2() && 
	   track->GetChiSquare()<MaxChisquareTr ){//MaXChisquare
	  
	  TrackCont.push_back(track);
	}    
	else{
	  delete track;
	}
      }
    }
  }
  
  for( int i=0; i<NumOfLayers/3; ++i ){
    for_each( CandContV[i].begin(), CandContV[i].end(), DeleteObject() );
    for_each( CandContX[i].begin(), CandContX[i].end(), DeleteObject() );
    for_each( CandContU[i].begin(), CandContU[i].end(), DeleteObject() );
  }

  {
    for_each( TrackContV.begin(), TrackContV.end(), DeleteObject() );
    for_each( TrackContX.begin(), TrackContX.end(), DeleteObject() );
    for_each( TrackContU.begin(), TrackContU.end(), DeleteObject() );
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }
  
  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp1() );    

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackCont[i];
	std::cout << std::setw(3) << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " DifVXU=" << track->GetDifVXU2()
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
#endif

#if 1
  // Delete Tracks about (Nhit1 = Nhit2  && chi1 < chi2)
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
      double chi=tp->GetChiSquare();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;
      double chi2=tp2->GetChiSquare();
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if((flag) && ((nh==nh2) && (chi<=chi2))){
        delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }    
#endif  

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": Before Deleted Sorting. #Tracks = " 
		<< nn << std::endl;
      
      //for( int i=0; i<100; ++i ){
      for( int i=0; i<nn; ++i ){
	TrLocalTrack *track=TrackCont[i];
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
  		TrackCont.end(), TrLTrackComp() );

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    if(nn>0){
      std::cout << funcname << ": After Sorting. #Tracks = " 
      		<< nn << std::endl;
      
      for( int i=0; i<nn; ++i ){
      //for( int i=0; i<100; ++i ){
	TrLocalTrack *track=TrackCont[i];
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
  
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nh=tp->GetNHit();
    for( int j=0; j<nh; ++j ) tp->GetHit(j)->setFlags();
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nh2=tp2->GetNHit(), flag=0;  
      for( int j=0; j<nh2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  // Clear Flags
  {
    int nbefore=TrackCont.size();
    for( int i=0; i<nbefore; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	tp->GetHit(j)->clearFlags();
      }
    }
  }  

  
  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nh=tp->GetNHit();
      for( int j=0; j<nh; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }
  
#if 0
  {
    int nn=TrackCont.size();
    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      if((track->GetChiSquare()>30)){
	std::cout << "" << std::endl;
	std::cout << std::setw(3) <<"ntrack = " << nn << " Track# =  " << i << " #Hits="
		  << std::setw(2) << track->GetNHit() 
		  << " ChiSqr=" << track->GetChiSquare()
		  << " DifVXU=" << track->GetDifVXU()
		  << std::endl;
	std::cout << "" << std::endl;

      }
    }
  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );


  return TrackCont.size();    
}


std::vector< std::vector<int> > 
makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 )
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
    index2=makeindex( ndim_org, minimumHit, ndim-1, index1+1 );
 
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
makeindex_below( int ndim_org, int maximumHit, int ndim, const int *index1 )
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
    index2=makeindex_below( ndim_org, maximumHit, ndim-1, index1+1 );
 
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
        if ((validHitNum <= maximumHit) && (validHitNum>0))
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      int size1=index.size();
    }
  }

  return index;
}


bool MakeHitCluster( const TrHitContainer & HC,
		     std::vector <TrHitCluster *> & Cont )
{
  int nh=HC.size(); 
  
  for( int i=0; i<nh; ++i ){
    TrHit *hit=HC[i];
    if( hit ){
      int multi = hit->GetPosSize();
      for (int m=0; m<multi; m++) {
	if( !(hit->rangecheck(m)) ) continue;
	double pos=hit->GetPos(m);	
	// 	double wp=hit->GetWirePosition();
	// 	double dl=hit->GetDriftLength(m);
	Cont.push_back( new TrHitCluster( new TrLTrackHit(hit,pos,m) ) );
      }
    }
  }

  return true;
}

TrLocalTrack *MakeTrack(  const std::vector < std::vector <TrHitCluster *> > &CandCont,
			  const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  TrLocalTrack *tp=new TrLocalTrack();

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    

  int n=CandCont.size();
  for( int i=0; i<n; ++i ){
    int m=combination[i];
    TrHitCluster *cluster=0;
    if(m>=0) cluster=CandCont[i][m];

#if 1
    if(cluster){
     int mm=cluster->NumberOfHits();
      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);
	  if(hitp) tp->AddHit( hitp );
      }
    }
#endif

#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif

#if 0
    /*****************************************************************/
    // Only SdcIn & BcOut Local Tracking !!                          
    /*****************************************************************/
    int mm = 0 ;
    //std::cout << "-------------------------------------" << std::endl;
    if(cluster){
     mm = cluster->NumberOfHits();
     if(mm>2) std::cout << "mm = " << mm << std::endl;

     double DL[mm] ;
     double DT[mm] ;
     int Layer[mm] ;

      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);

	DL[j] = hitp->GetDriftLength();
	DT[j] = hitp->GetDriftTime();
	Layer[j] = hitp->GetLayer();

	if(mm<2 || !((Layer[j]>0 && Layer[j]<11) || (Layer[j]>112 && Layer[j]<125))){
	  if((Layer[j]>=1 && Layer[j]<=4) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=5 && Layer[j]<=10) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=113 && Layer[j]<=118) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=119 && Layer[j]<=124) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	}
      }
      //std::cout << "***********************************" << std::endl;
      if(mm==2 && ((Layer[0]>0 && Layer[0]<11) || (Layer[0]>112 && Layer[0]<125))){
	int Ok = 0 ;
	/*	
	std::cout << "DL[ " << Layer[0] << "] = " << DL[0] 
		  << " DL[ " << Layer[1] << "] = " << DL[1] << std::endl; 
	*/
	//std::cout << "***********************************" << std::endl;
	if(Layer[0] ==1 || Layer[0]==3 || Layer[0]==113 || Layer[0]==115 || Layer[0]==117){
	  if(((DL[0]+DL[1])>0.7) && ((DL[0]+DL[1])<2.3)) Ok = 1 ;
	}  
	if(Layer[0] ==5 || Layer[0]==7 || Layer[0]==9 || Layer[0]==119 || Layer[0]==121 || Layer[0]==123){
	  if(((DL[0]+DL[1])>1.5) && ((DL[0]+DL[1])<3.5)) Ok = 1 ;
	}
	if(Ok){
	  for(int j=0; j<mm; ++j ){
	    TrLTrackHit *hitp=cluster->GetHit(j);
	    if(hitp) tp->AddHit( hitp );
	  }
	}
      }
    }
#endif

  }
  return tp;
}

