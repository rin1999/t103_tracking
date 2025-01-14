/*
  TrAnalyzer.hh

  2024/11  K.Shirotori
*/

#ifndef TrAnalyzer_h 
#define TrAnalyzer_h 1

#include "DetectorInfo.hh"
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
  TrHitContainer BFTHC[NumOfLayersBFT+1];

   // TrHit Cluserized 
  TrHitContainer BFTCHC[NumOfLayersBFT+1];

  //TrCluster
  TrClusterContainer BFTClCont;

  // TrLocalTrack
  std::vector <TrLocalTrack *> TrackBFTCol;

  // hit layers of local track
  //std::vector <std::vector<int>*> hitlayersCont;

public:
  bool DecodeBFTRawHits( RawData *rawData );

  inline const TrHitContainer & GetBFTHC( int layer ) const;
  // TrHit Cluserized 
  inline const TrHitContainer & GetBFTCHC( int layer ) const;

  int GetNClustersBFT( void ) const { return BFTClCont.size(); };
  inline TrCluster * GetClusterBFT( int i ) const;

  bool TrackSearchBFT( void );

  int GetNtracksBFT( void ) const  { return TrackBFTCol.size(); }

  inline TrLocalTrack * GetTrackBFT( int i ) const;

  //std::vector<std::vector<int>*> GetHitLayers() const { return hitlayersCont; }

  bool ReCalcBFTClusters( bool applyRecursively=false );

  bool ReCalcBFTHits( bool applyRecursively=false ); 

  bool ReCalcTrackBFT( bool applyRecursively=false ); 

  bool ReCalcAll( void );

private:
  void clearBFTHits( void );
  void clearBFTCHits( void );
  void clearTracksBFT( void );

  static int MakeUpClusters( const TrHitContainer & HitCont,
			     TrClusterContainer & ClusterCont,
			     double maxTimeDif );

public:
  void resetTracksTr( void ) { clearTracksBFT(); }
};

inline const TrHitContainer & TrAnalyzer::GetBFTHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBFT ) layer=0;
  return BFTHC[layer];
}

inline const TrHitContainer & TrAnalyzer::GetBFTCHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBFT ) layer=0;
  return BFTCHC[layer];
}

inline TrCluster * TrAnalyzer::GetClusterBFT( int i ) const
{
  if( i>=0 && i<BFTClCont.size() )
    return BFTClCont[i];
  else
    return 0;
}

inline TrLocalTrack * TrAnalyzer::GetTrackBFT( int i ) const
{
  if( i>=0 && i<TrackBFTCol.size() )
    return TrackBFTCol[i];
  else
    return 0;
}

#endif 
