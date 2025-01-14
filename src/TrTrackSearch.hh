/*
  TrTrackSearch.hh

  2019/2  K.Shirotori
*/

#ifndef TrTrackSearch_h
#define TrTrackSearch_h 1

#include "TrAnalyzer.hh"

#include <vector>

struct TrPairPlaneInfo;
class TrPairHitCluster;
class TrLocalTrack;
class TrLTrackHit;

int LocalTrackSearch( const TrHitContainer * HC,
		      std::vector <TrLocalTrack *> &TrackCont,
		      int NumOfLayers, int MinNumOfHits);

bool MakePairPlaneHitCluster( const TrHitContainer & HC1,
			      const TrHitContainer & HC2,
			      double CellSize,
			      std::vector <TrPairHitCluster *> & Cont );


bool MakeUnPairPlaneHitCluster( const TrHitContainer & HC,
				std::vector <TrPairHitCluster *> & Cont );


TrLocalTrack *MakeTrack( const std::vector < std::vector <TrPairHitCluster *> > &CandCont,
			 const int *combination );

std::vector< std::vector<int> > makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 ); 

#endif
