/*
  DCTrackSearch.hh
*/

#ifndef DCTrackSearch_h
#define DCTrackSearch_h 1

#include "DCAnalyzer.hh"

#include <vector>

struct DCPairPlaneInfo;
class DCPairHitCluster;
class DCLocalTrack;
class DCLTrackHit;
class MWPCCluster;

//For BC1&2
int MWPCLocalTrackSearch( const DCHitContainer * HC,
			  std::vector <DCLocalTrack *> & TrackCont);

//For BC3&4 SDC1&2
int LocalTrackSearch( const DCHitContainer * HC,  
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> & TrackCont,
		      int MinNumOfHits=6 );

//For SDC3&4
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont );


int LocalTrackSearchBcOutSdcIn( const DCHitContainer * BcHC,  
                                 const DCPairPlaneInfo * BcPpInfo,
                                 const DCHitContainer * SdcHC,  
                                 const DCPairPlaneInfo * SdcPpInfo,
                                 int BcNpp, int SdcNpp,
                                 std::vector <DCLocalTrack *> &TrackCont,
                                 int MinNumOfHits=18 );

bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont );

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont );

bool MakeMWPCPairPlaneHitCluster( const DCHitContainer & HC,  
				  std::vector <DCPairHitCluster *> & Cont );

DCLocalTrack *MakeTrack( const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
			 const int *combination );

std::vector< std::vector<int> > makeindex( int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_BcIn( int ndim_org, int minimumHit, int ndim, const int *index1 ); 

// bool MakeSdc4HitCluster( const DCHitContainer & HC1,
// 			 const DCHitContainer & HC2,
// 			 const DCHitContainer & HC3,
// 			 const DCHitContainer & HC4,
// 			 std::vector <DC4HitCluster *> & Cont );

// bool MakeSdc4HitCluster( const DCHitContainer & HC1,
// 			 const DCHitContainer & HC2,
// 			 const DCHitContainer & HC3,
// 			 const DCHitContainer & HC4,
// 			 const DCHitContainer & HC5,
// 			 const DCHitContainer & HC6,
// 			 std::vector <DC4HitCluster *> & Cont );

// DCLocalTrack *MakeTrackAndCheck( DCPairHitCluster * cluster1, DC4HitCluster * cluster2,
// 				 DCLTrackHit * hit1, DCLTrackHit * hit2 );

// DCLocalTrack *CombineTrackXY( DCLocalTrack *trX, DCLocalTrack *trY );

#endif
