/*
  TrAnalyzer.hh

  2019/2  K.Shirotori
*/

#ifndef TrAnalyzer_h 
#define TrAnalyzer_h 1

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

class RawData;
class TrHit;
class TrLocalTrack;
class TrCluster;

typedef std::vector <TrHit *> TrHitContainer;
typedef std::vector <TrCluster *> TrClusterContainer;

class TrAnalyzer
{
public:
  TrAnalyzer();
  ~TrAnalyzer();
private:
  TrAnalyzer( const TrAnalyzer & );
  TrAnalyzer & operator = ( const TrAnalyzer & );

private:
  // TrHit
  TrHitContainer TrHC[NumOfLayersSFT+1];
  // TrHit Cluserized 
  TrHitContainer TrCHC[NumOfLayersSFT+1];

  //TrCluster
  TrClusterContainer TrClCont;

  // TrLocalTrack
  std::vector <TrLocalTrack *> TrackTrCol;

public:
  bool DecodeRawHits( RawData *rawData );

  inline const TrHitContainer & GetTrHC( int layer ) const;
  // TrHit Cluserized 
  inline const TrHitContainer & GetTrCHC( int layer ) const;

  int GetNClustersTr( void ) const { return TrClCont.size(); };
  inline TrCluster * GetClusterTr( int i ) const;

  bool TrackSearchTr( void );

  int GetNtracksTr( void ) const  { return TrackTrCol.size(); }

  inline TrLocalTrack * GetTrackTr( int i ) const;

  bool ReCalcTrClusters( bool applyRecursively=false );

  bool ReCalcTrHits( bool applyRecursively=false ); 

  bool ReCalcTrackTr( bool applyRecursively=false ); 

  bool ReCalcAll( void );

private:
  void clearTrHits( void );
  void clearTracksTr( void );

  static int MakeUpClusters( const TrHitContainer & HitCont,
			     TrClusterContainer & ClusterCont,
			     double maxTimeDif );

public:
  void resetTracksTr( void ) { clearTracksTr(); }
};

inline const TrHitContainer & TrAnalyzer::GetTrHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=0;
  return TrHC[layer];
}

inline const TrHitContainer & TrAnalyzer::GetTrCHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=0;
  return TrCHC[layer];
}

inline TrCluster * TrAnalyzer::GetClusterTr( int i ) const
{
  if( i>=0 && i<TrClCont.size() )
    return TrClCont[i];
  else
    return 0;
}

inline TrLocalTrack * TrAnalyzer::GetTrackTr( int i ) const
{
  if( i>=0 && i<TrackTrCol.size() )
    return TrackTrCol[i];
  else
    return 0;
}

#endif 
