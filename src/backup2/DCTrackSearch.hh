/*
  DCTrackSearch.hh

  2018/12  K.Shirotori
*/

#ifndef DCTrackSearch_h
#define DCTrackSearch_h 1

#include "DCAnalyzer.hh"

#include <vector>

struct DCPairPlaneInfo;
class DCPairHitCluster;
class DCLocalTrack;
class DCLTrackHit;

int LocalTrackSearch( const DCHitContainer * HC,
		      std::vector <DCLocalTrack *> &TrackCont,
		      int NumOfLayers, int MinNumOfHits);

bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont );


bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont );


DCLocalTrack *MakeTrack( const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
			 const int *combination );

std::vector< std::vector<int> > makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 ); 

#endif
