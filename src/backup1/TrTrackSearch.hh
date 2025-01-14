/*
  TrTrackSearch.hh

  2012/5  K.Shirotori
*/

#ifndef TrTrackSearch_h
#define TrTrackSearch_h 1

#include "TrAnalyzer.hh"

#include <vector>

class TrHitCluster;
class TrLocalTrack;
class TrLTrackHit;

//Full combination by Linear fitting
int LocalTrackSearch( const TrHitContainer * HC,
		      std::vector <TrLocalTrack *> &TrackCont,
		      int NumOfLayers, int MinNumOfHits,
		      bool FlgMinSegs = false);

//Full combination by Quadra fitting
int LocalTrackSearchQ( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont,
		       int NumOfLayers, int MinNumOfHits);

//Full combination by Hyperbola fitting
int LocalTrackSearchH( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont,
		       int NumOfLayers, int MinNumOfHits);

//VXU tracking by Linear fitting
int LocalTrackSearch2( const TrHitContainer * HC,
		       std::vector <TrLocalTrack *> &TrackCont, 
		       int NumOfLayers, int MinNumOfHits);

//VXU tracking by by Quadra -> Quadra fitting
int LocalTrackSearchQ2( const TrHitContainer * HC,
			std::vector <TrLocalTrack *> &TrackCont,
			int NumOfLayers, int MinNumOfHits);

// //VXU tracking by by  Linear -> Quadra fitting
// int LocalTrackSearchQ2( const TrHitContainer * HC,
// 			std::vector <TrLocalTrack *> &TrackCont,
// 			int NumOfLayers, int MinNumOfHits);

std::vector< std::vector<int> > makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_below( int ndim_org, int maximumHit, int ndim, const int *index1 ); 

bool MakeHitCluster( const TrHitContainer & HC,  
		     std::vector <TrHitCluster *> & Cont );

TrLocalTrack *MakeTrack( const std::vector < std::vector <TrHitCluster *> > &CandCont,
			 const int *combination );

#endif
