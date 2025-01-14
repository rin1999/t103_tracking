/*
  DCAnalyzer.hh
*/

#ifndef DCAnalyzer_h 
#define DCAnalyzer_h 1

#include "DetectorID.hh"
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
  DCHitContainer TestDCHC[NumOfLayersDC+1];

  // DCLocalTrack
  std::vector <DCLocalTrack *> TrackTestDCCol;

public:
  bool DecodeRawHits( RawData *rawData );

  inline const DCHitContainer & GetTestDCHC( int layer ) const;

  bool TrackSearchTestDC( void );

  int GetNtracksTestDC( void ) const  { return TrackTestDCCol.size(); }

  inline DCLocalTrack * GetTrackTestDC( int i ) const;

  bool ReCalcTestDCHits( bool applyRecursively=false ); 

  bool ReCalcTrackTestDC( bool applyRecursively=false ); 

  bool ReCalcAll( void );

private:
  void clearTestDCHits( void );
  void clearTracksTestDC( void );

public:
  void resetTracksTestDC( void ) { clearTracksTestDC(); }
};

inline const DCHitContainer & DCAnalyzer::GetTestDCHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersDC ) layer=0;
  return TestDCHC[layer];
}

inline DCLocalTrack * DCAnalyzer::GetTrackTestDC( int i ) const
{
  if( i>=0 && i<TrackTestDCCol.size() )
    return TrackTestDCCol[i];
  else
    return 0;
}

#endif 
