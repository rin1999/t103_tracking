/*
  DCAnalyzer.hh
*/

#ifndef DCAnalyzer_h 
#define DCAnalyzer_h 1

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class DCHit;
class DCLocalTrack;
class K18Track;
class SksTrack;
class RawData;
class MWPCCluster;
//class SimuData;

class Hodo1Hit;
class Hodo2Hit;

typedef std::vector <DCHit *> DCHitContainer;
typedef std::vector <MWPCCluster*> MWPCClusterContainer;

typedef std::vector <Hodo1Hit *> Hodo1HitContainer;
typedef std::vector <Hodo2Hit *> Hodo2HitContainer;

class DCAnalyzer
{
public:
  DCAnalyzer();
  ~DCAnalyzer();
private:
  DCAnalyzer( const DCAnalyzer & );
  DCAnalyzer & operator = ( const DCAnalyzer & );

private:
  // DCHit
  DCHitContainer TempBcInHC[NumOfLayersBcIn+1];
  DCHitContainer BcInHC[NumOfLayersBcIn+1], BcOutHC[NumOfLayersBcOut+1];
  DCHitContainer SdcInHC[NumOfLayersSdcIn+1], SdcOutHC[NumOfLayersSdcOut+1];
  /*New Added by K.Miwa*/
  DCHitContainer VtxPoint;

  //MWPC clustering
  MWPCClusterContainer MWPCClCont[NumOfLayersBcIn+1];

  // DCLocalTrack
  std::vector <DCLocalTrack *> TrackBcInCol, TrackBcOutCol;
  std::vector <DCLocalTrack *> TrackSdcInCol, TrackSdcOutCol;

  //DCLocalTrack -> for SdcOut Tracking
  std::vector <DCLocalTrack *> TrackSdcOutCol1, TrackSdcOutCol2;

  // K18Track
  std::vector <K18Track *> K18TrackCol;

  // SksTrack
  std::vector <SksTrack *> SksTrackCol;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  bool DecodeRawHits( RawData *rawData );
  /* New Added by K.Miwa */
  //bool DecodeSimuHits( SimuData *simuData );

  inline const DCHitContainer & GetTempBcInHC( int layer ) const;
  inline const DCHitContainer & GetBcInHC( int layer ) const;
  inline const DCHitContainer & GetBcOutHC( int layer ) const;
  inline const DCHitContainer & GetSdcInHC( int layer ) const;
  inline const DCHitContainer & GetSdcOutHC( int layer ) const;

  bool TrackSearchBcIn( void );
  bool TrackSearchBcIn( const std::vector<std::vector<DCHitContainer> >& hc );
  bool TrackSearchBcOut( void );
  bool TrackSearchSdcIn( void );
  bool TrackSearchSdcOut( void );

  int GetNtracksBcIn( void ) const  { return TrackBcInCol.size(); }
  int GetNtracksBcOut( void ) const { return TrackBcOutCol.size(); }
  int GetNtracksSdcIn( void ) const  { return TrackSdcInCol.size(); }
  int GetNtracksSdcOut( void ) const { return TrackSdcOutCol.size(); }

  inline DCLocalTrack * GetTrackBcIn( int i ) const;
  inline DCLocalTrack * GetTrackBcOut( int i ) const;
  inline DCLocalTrack * GetTrackSdcIn( int i ) const;
  inline DCLocalTrack * GetTrackSdcOut( int i ) const;

  bool TrackSearchK18( void );
  bool TrackSearchSks( void );

  int GetNTracksK18( void ) const { return K18TrackCol.size(); }
  int GetNTracksSks( void ) const { return SksTrackCol.size(); }

  inline K18Track  * GetK18Track( int i ) const;
  inline SksTrack * GetSksTrack( int i ) const;

  int GetNClustersMWPC( int layer ) const { return MWPCClCont[layer].size(); };
  inline const MWPCClusterContainer & GetClusterMWPC( int layer ) const;

  bool ReCalcDCHits( bool applyRecursively=false ); 

  bool ReCalcTrackBcIn( bool applyRecursively=false ); 
  bool ReCalcTrackBcOut( bool applyRecursively=false ); 
  bool ReCalcTrackSdcIn( bool applyRecursively=false ); 
  bool ReCalcTrackSdcOut( bool applyRecursively=false ); 

  bool ReCalcK18Track( bool applyRecursively=false ); 
  bool ReCalcSksTrack( bool applyRecursively=false ); 

  bool ReCalcAll( void );

  // Optional Extension by Miwa
private:
  std::vector <DCLocalTrack *> TrackBcOutSdcInCol;
public:
  bool TrackSearchBcOutSdcIn( void );
  int GetNtracksBcOutSdcIn( void ) const { return TrackBcOutSdcInCol.size(); }
  inline DCLocalTrack * GetTrackBcOutSdcIn( int i ) const;
private:
  void clearTracksBcOutSdcIn( void ) ;
  // End of Extension

private:
  void clearDCHits( void );
  void clearVtxHits( void );
  void clearTracksBcIn( void );
  void clearTracksBcOut( void );
  void clearTracksSdcIn( void );
  void clearTracksSdcOut( void );
  void clearK18Tracks( void );
  void clearSksTracks( void );

  static int MakeUpMWPCClusters( const DCHitContainer & HitCont,
				 MWPCClusterContainer & ClusterCont,
				 double maxTimeDif );
public:
  void resetTracksBcIn( void ) { clearTracksBcIn(); }
  void resetTracksBcOut( void ) { clearTracksBcOut(); }
  void resetTracksSdcIn( void ) { clearTracksSdcIn(); }
  void resetTracksSdcOut( void ) { clearTracksSdcOut(); }
  void resetTracksBcOutSdcIn( void ) { clearTracksBcOutSdcIn(); }

  // for DS
// public:
//   std::size_t DSSize( void ) const;
//   std::size_t DSSave( unsigned int *bufp );
//   bool DSRestore( unsigned int *bufp );

};


inline const DCHitContainer & DCAnalyzer::GetTempBcInHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return TempBcInHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetBcInHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return BcInHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetBcOutHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcOut ) layer=0;
  return BcOutHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetSdcInHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcIn ) layer=0;
  return SdcInHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetSdcOutHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcOut ) layer=0;
  return SdcOutHC[layer];
}

inline DCLocalTrack * DCAnalyzer::GetTrackBcIn( int i ) const
{
  if( i>=0 && i<TrackBcInCol.size() )
    return TrackBcInCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackBcOut( int i ) const
{
  if( i>=0 && i<TrackBcOutCol.size() )
    return TrackBcOutCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackSdcIn( int i ) const
{
  if( i>=0 && i<TrackSdcInCol.size() )
    return TrackSdcInCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackSdcOut( int i ) const
{
  if( i>=0 && i<TrackSdcOutCol.size() )
    return TrackSdcOutCol[i];
  else
    return 0;
}

inline K18Track * DCAnalyzer::GetK18Track( int i ) const
{
  if( i>=0 && i<K18TrackCol.size() )
    return K18TrackCol[i];
  else
    return 0;
}

inline SksTrack * DCAnalyzer::GetSksTrack( int i ) const
{
  if( i>=0 && i<SksTrackCol.size() )
    return SksTrackCol[i];
  else
    return 0;
}

// Optional Extension by Miwa
inline DCLocalTrack * DCAnalyzer::GetTrackBcOutSdcIn( int i ) const
{
  if( i>=0 && i<TrackBcOutSdcInCol.size() )
    return TrackBcOutSdcInCol[i];
  else
    return 0;
}
// End of Extension

inline const MWPCClusterContainer & DCAnalyzer::GetClusterMWPC( int layer ) const
{ 
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return MWPCClCont[layer];
}

#endif 
