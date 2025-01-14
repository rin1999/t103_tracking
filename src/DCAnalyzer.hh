/*
  DCAnalyzer.hh
*/

#ifndef DCAnalyzer_h 
#define DCAnalyzer_h 1

#include "DetectorInfo.hh"
#include "ThreeVector.hh"
#include <vector>

class DCHit;
class DCLocalTrack;
class RawData;

typedef std::vector <DCHit *> DCHitContainer;

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
   DCHitContainer BDCHC[NumOfLayersBDC+1];
   DCHitContainer KLDCHC[NumOfLayersKLDC+1];
   
   // DCLocalTrack
   std::vector <DCLocalTrack *> TrackBDCCol;
   std::vector <DCLocalTrack *> TrackKLDCCol;
   
public:
   bool DecodeBDCRawHits( RawData *rawData );
   bool DecodeKLDCRawHits( RawData *rawData );
   
   inline const DCHitContainer & GetBDCHC( int layer ) const;
   inline const DCHitContainer & GetKLDCHC( int layer ) const;
   
   bool TrackSearchBDC( void );
   bool TrackSearchKLDC( void );
   
   int GetNtracksBDC( void ) const  { return TrackBDCCol.size(); }
   int GetNtracksKLDC( void ) const  { return TrackKLDCCol.size(); }
   
   inline DCLocalTrack * GetTrackBDC( int i ) const;
   inline DCLocalTrack * GetTrackKLDC( int i ) const;
   
   bool ReCalcBDCHits( bool applyRecursively=false ); 
   bool ReCalcKLDCHits( bool applyRecursively=false ); 
   
   bool ReCalcTrackBDC( bool applyRecursively=false ); 
   bool ReCalcTrackKLDC( bool applyRecursively=false ); 
   
   bool ReCalcAll( void );

private:
  void clearBDCHits( void );
  void clearKLDCHits( void );
  void clearTracksBDC( void );
  void clearTracksKLDC( void );

public:
  void resetTracksBDC( void ) { clearTracksBDC(); }
  void resetTracksKLDC( void ) { clearTracksKLDC(); }
};

inline const DCHitContainer & DCAnalyzer::GetBDCHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBDC ) layer=0;
  return BDCHC[layer];
}

inline const DCHitContainer & DCAnalyzer::GetKLDCHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersKLDC ) layer=0;
  return KLDCHC[layer];
}

inline DCLocalTrack * DCAnalyzer::GetTrackBDC( int i ) const
{
  if( i>=0 && i<TrackBDCCol.size() )
    return TrackBDCCol[i];
  else
    return 0;
}

inline DCLocalTrack * DCAnalyzer::GetTrackKLDC( int i ) const
{
  if( i>=0 && i<TrackKLDCCol.size() )
    return TrackKLDCCol[i];
  else
    return 0;
}

#endif 
