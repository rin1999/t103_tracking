/*
  TrAnalyzer.hh

  2016/2  K.Shirotori
*/

#ifndef TrAnalyzer_h 
#define TrAnalyzer_h 1

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

class TrHit;
class TrLocalTrack;

// class BeamTrack;
// class PreInTrack;
// class PreOutTrack;
// class PreOut2Track;
// class Scat1Track;
// class Scat2Track;
// class Scat2ATrack;
// class Scat2BTrack;
// class Scat3Track;

class RawData;
class s_BeamRawHit;
class s_ScatRawHit;

typedef std::vector <TrHit *> TrHitContainer;

class TrAnalyzer
{
public:
  TrAnalyzer();
  ~TrAnalyzer();
private:
  TrAnalyzer( const TrAnalyzer & );
  TrAnalyzer & operator = ( const TrAnalyzer & );

private:
  //Simple tracking

  TrHitContainer bSSDTHC[NumOfLayersbSSD+1];
  TrHitContainer sSSD1THC[NumOfLayerssSSD1+1];
  TrHitContainer sSSD2THC[NumOfLayerssSSD2+1];
  // TrHitContainer BFTTHC[NumOfLayersBeamT+1];
  // TrHitContainer SFTTHC[NumOfLayersSFT+1];
  // TrHitContainer IT1THC[NumOfLayersIT1+1];
  // TrHitContainer IT2RTHC[PlMaxIT2R+1];
  // TrHitContainer IT2LTHC[PlMaxIT2L+1];
  // TrHitContainer ST1THC[NumOfLayersST1+1];
  // TrHitContainer ST2THC[NumOfLayersST2+1];
  // TrHitContainer ScatInTHC[NumOfLayersScatInT+1];

  std::vector <TrLocalTrack *> TrackbSSDTCol;
  std::vector <TrLocalTrack *> TracksSSD2TCol;
  std::vector <TrLocalTrack *> TracksSSD1TCol;
  // std::vector <TrLocalTrack *> TrackBFTTCol;
  // std::vector <TrLocalTrack *> TrackSFTTCol;
  // std::vector <TrLocalTrack *> TrackIT1TCol;
  // std::vector <TrLocalTrack *> TrackIT2RTCol;
  // std::vector <TrLocalTrack *> TrackIT2LTCol;
  // std::vector <TrLocalTrack *> TrackST1TCol;
  // std::vector <TrLocalTrack *> TrackST2TCol;
  // std::vector <TrLocalTrack *> TrackScatInTCol;

  // std::vector <BeamTrack *> BeamTrackCol;
  // std::vector <PreInTrack *> PreInTrackCol;
  // std::vector <PreOutTrack *> PreOutTrackCol;
  // std::vector <PreOut2Track *> PreOut2TrackCol;
  // std::vector <Scat1Track *> Scat1TrackCol;
  // std::vector <Scat2Track *> Scat2TrackCol;
  // std::vector <Scat2ATrack *> Scat2ATrackCol;
  // std::vector <Scat2BTrack *> Scat2BTrackCol;
  // std::vector <Scat3Track *> Scat3TrackCol;

public:
  bool DecodeRawHits( RawData *rawData );

  bool DecodesBRawHits( s_BeamRawHit *sbeamRaw, int type );
  bool DecodesSRawHits( s_ScatRawHit *sscatRaw, int type );

  inline const TrHitContainer & GetbSSDTHC( int layer ) const;
  inline const TrHitContainer & GetsSSD1THC( int layer ) const;
  inline const TrHitContainer & GetsSSD2THC( int layer ) const;
  // inline const TrHitContainer & GetBFTTHC( int layer ) const;
  // inline const TrHitContainer & GetSFTTHC( int layer ) const;
  // inline const TrHitContainer & GetIT1THC( int layer ) const;
  // inline const TrHitContainer & GetIT2RTHC( int layer ) const;
  // inline const TrHitContainer & GetIT2LTHC( int layer ) const;
  // inline const TrHitContainer & GetST1THC( int layer ) const;
  // inline const TrHitContainer & GetST2THC( int layer ) const;
  // inline const TrHitContainer & GetScatInTHC( int layer ) const;

  bool TrackSearchbSSDT( void );
  bool TrackSearchsSSD1T( void );
  bool TrackSearchsSSD2T( void );
  // bool TrackSearchBFTT( void );
  // bool TrackSearchSFTT( void );
  // bool TrackSearchIT1T( void );
  // bool TrackSearchIT2RT( void );
  // bool TrackSearchIT2LT( void );
  // bool TrackSearchST1T( void );
  // bool TrackSearchST2T( void );
  // bool TrackSearchScatInT( void );

  int GetNtracksbSSDT( void ) const  { return TrackbSSDTCol.size(); }
  int GetNtrackssSSD1T( void ) const  { return TracksSSD1TCol.size(); }
  int GetNtrackssSSD2T( void ) const  { return TracksSSD2TCol.size(); }
  // int GetNtracksBFTT( void ) const  { return TrackBFTTCol.size(); }
  // int GetNtracksSFTT( void ) const  { return TrackSFTTCol.size(); }
  // int GetNtracksIT1T( void ) const  { return TrackIT1TCol.size(); }
  // int GetNtracksIT2RT( void ) const  { return TrackIT2RTCol.size(); }
  // int GetNtracksIT2LT( void ) const  { return TrackIT2LTCol.size(); }
  // int GetNtracksST1T( void ) const  { return TrackST1TCol.size(); }
  // int GetNtracksST2T( void ) const  { return TrackST2TCol.size(); }
  // int GetNtracksScatInT( void ) const  { return TrackScatInTCol.size(); }

  inline TrLocalTrack * GetTrackbSSDT( int i ) const;
  inline TrLocalTrack * GetTracksSSD1T( int i ) const;
  inline TrLocalTrack * GetTracksSSD2T( int i ) const;
  // inline TrLocalTrack * GetTrackBFTT( int i ) const;
  // inline TrLocalTrack * GetTrackSFTT( int i ) const;
  // inline TrLocalTrack * GetTrackIT1T( int i ) const;
  // inline TrLocalTrack * GetTrackIT2RT( int i ) const;
  // inline TrLocalTrack * GetTrackIT2LT( int i ) const;
  // inline TrLocalTrack * GetTrackST1T( int i ) const;
  // inline TrLocalTrack * GetTrackST2T( int i ) const;
  // inline TrLocalTrack * GetTrackScatInT( int i ) const;

  // bool TrackSearchBeam( ThreeVector vertex );
  // bool TrackSearchPreIn( double IniP );
  // bool TrackSearchPreOut( double IniP );
  // bool TrackSearchPreOut2( double IniP );
  // bool TrackSearchScat1( double IniP, ThreeVector vertex );
  // bool TrackSearchScat2( double IniP, ThreeVector vertex, s_ScatRawHit *sscatRaw );
  // bool TrackSearchScat3( double IniP, ThreeVector vertex );

  // int GetNTracksBeam( void ) const { return BeamTrackCol.size(); }
  // int GetNTracksPreIn( void ) const { return PreInTrackCol.size(); }
  // int GetNTracksPreOut( void ) const { return PreOutTrackCol.size(); }
  // int GetNTracksPreOut2( void ) const { return PreOut2TrackCol.size(); }
  // int GetNTracksScat1( void ) const { return Scat1TrackCol.size(); }
  // int GetNTracksScat2( void ) const { return Scat2TrackCol.size(); }
  // int GetNTracksScat2A( void ) const { return Scat2ATrackCol.size(); }
  // int GetNTracksScat2B( void ) const { return Scat2BTrackCol.size(); }
  // int GetNTracksScat3( void ) const { return Scat3TrackCol.size(); }

  // inline BeamTrack * GetBeamTrack( int i ) const;
  // inline PreInTrack * GetPreInTrack( int i ) const;
  // inline PreOutTrack * GetPreOutTrack( int i ) const;
  // inline PreOut2Track * GetPreOut2Track( int i ) const;
  // inline Scat1Track * GetScat1Track( int i ) const;
  // inline Scat2Track * GetScat2Track( int i ) const;
  // inline Scat2ATrack * GetScat2ATrack( int i ) const;
  // inline Scat2BTrack * GetScat2BTrack( int i ) const;
  // inline Scat3Track * GetScat3Track( int i ) const;

  // bool ReCalcTrHits( bool applyRecursively=false ); 
  // bool ReCalcTrackBFTT( bool applyRecursively=false ); 
  // bool ReCalcTrackSFTT( bool applyRecursively=false ); 
  // bool ReCalcTrackIT1T( bool applyRecursively=false ); 
  // bool ReCalcTrackIT2RT( bool applyRecursively=false ); 
  // bool ReCalcTrackIT2LT( bool applyRecursively=false ); 
  // bool ReCalcTrackST1T( bool applyRecursively=false ); 
  // bool ReCalcTrackST2T( bool applyRecursively=false ); 
  // bool ReCalcTrackScatInT( bool applyRecursively=false ); 

  // bool ReCalcBeamTrack( bool applyRecursively=false ); 
  // bool ReCalcPreInTrack( bool applyRecursively=false ); 
  // bool ReCalcPreOutTrack( bool applyRecursively=false ); 
  // bool ReCalcPreOut2Track( bool applyRecursively=false ); 
  // bool ReCalcScat1Track( bool applyRecursively=false ); 
  // bool ReCalcScat2Track( bool applyRecursively=false ); 
  // bool ReCalcScat2ATrack( bool applyRecursively=false ); 
  // bool ReCalcScat2BTrack( bool applyRecursively=false ); 
  // bool ReCalcScat3Track( bool applyRecursively=false ); 

  // bool ReCalcAll( void );

private:
  void clearTrHits( void );
  void clearTracksbSSDT( void );
  void clearTrackssSSD1T( void );
  void clearTrackssSSD2T( void );
  // void clearTracksBFTT( void );
  // void clearTracksSFTT( void );
  // void clearTracksIT1T( void );
  // void clearTracksIT2RT( void );
  // void clearTracksIT2LT( void );
  // void clearTracksST1T( void );
  // void clearTracksST2T( void );
  // void clearTracksScatInT( void );
  // void clearBeamTracks( void );
  // void clearPreInTracks( void );
  // void clearPreOutTracks( void );
  // void clearPreOut2Tracks( void );
  // void clearScat1Tracks( void );
  // void clearScat2Tracks( void );
  // void clearScat2ATracks( void );
  // void clearScat2BTracks( void );
  // void clearScat3Tracks( void );

public:
  void resetTracksbSSDT( void ) { clearTracksbSSDT(); }
  void resetTrackssSSD1T( void ) { clearTrackssSSD1T(); }
  void resetTrackssSSD2T( void ) { clearTrackssSSD2T(); }
  // void resetTracksBFTT( void ) { clearTracksBFTT(); }
  // void resetTracksSFTT( void ) { clearTracksSFTT(); }
  // void resetTracksIT1T( void ) { clearTracksIT1T(); }
  // void resetTracksIT2RT( void ) { clearTracksIT2RT(); }
  // void resetTracksIT2LT( void ) { clearTracksIT2LT(); }
  // void resetTracksST1T( void ) { clearTracksST1T(); }
  // void resetTracksST2T( void ) { clearTracksST2T(); }
  // void resetTracksScatInT( void ) { clearTracksScatInT(); }
};

inline const TrHitContainer & TrAnalyzer::GetbSSDTHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersbSSD ) layer=0;
  return bSSDTHC[layer];
}

inline const TrHitContainer & TrAnalyzer::GetsSSD1THC( int layer ) const
{
  if( layer<0 || layer>NumOfLayerssSSD1 ) layer=0;
  return sSSD1THC[layer];
}

inline const TrHitContainer & TrAnalyzer::GetsSSD2THC( int layer ) const
{
  if( layer<0 || layer>NumOfLayerssSSD2 ) layer=0;
  return sSSD2THC[layer];
}

// inline const TrHitContainer & TrAnalyzer::GetBFTTHC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersBFT ) layer=0;
//   return BFTTHC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetSFTTHC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersSFT ) layer=0;
//   return SFTTHC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetIT1THC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersIT1 ) layer=0;
//   return IT1THC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetIT2RTHC( int layer ) const
// {
//   if( layer<0 || layer>PlMaxIT2R ) layer=0;
//   return IT2RTHC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetIT2LTHC( int layer ) const
// {
//   if( layer<0 || layer>PlMaxIT2L ) layer=0;
//   return IT2LTHC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetST1THC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersST1 ) layer=0;
//   return ST1THC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetST2THC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersST1 ) layer=0;
//   return ST2THC[layer];
// }

// inline const TrHitContainer & TrAnalyzer::GetScatInTHC( int layer ) const
// {
//   if( layer<0 || layer>NumOfLayersScatInT ) layer=0;
//   return ScatInTHC[layer];
// }

inline TrLocalTrack * TrAnalyzer::GetTrackbSSDT( int i ) const
{
  if( i>=0 && i<TrackbSSDTCol.size() )
    return TrackbSSDTCol[i];
  else
    return 0;
}

inline TrLocalTrack * TrAnalyzer::GetTracksSSD1T( int i ) const
{
  if( i>=0 && i<TracksSSD1TCol.size() )
    return TracksSSD1TCol[i];
  else
    return 0;
}

inline TrLocalTrack * TrAnalyzer::GetTracksSSD2T( int i ) const
{
  if( i>=0 && i<TracksSSD2TCol.size() )
    return TracksSSD2TCol[i];
  else
    return 0;
}

// inline TrLocalTrack * TrAnalyzer::GetTrackBFTT( int i ) const
// {
//   if( i>=0 && i<TrackBFTTCol.size() )
//     return TrackBFTTCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackSFTT( int i ) const
// {
//   if( i>=0 && i<TrackSFTTCol.size() )
//     return TrackSFTTCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackIT1T( int i ) const
// {
//   if( i>=0 && i<TrackIT1TCol.size() )
//     return TrackIT1TCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackIT2RT( int i ) const
// {
//   if( i>=0 && i<TrackIT2RTCol.size() )
//     return TrackIT2RTCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackIT2LT( int i ) const
// {
//   if( i>=0 && i<TrackIT2LTCol.size() )
//     return TrackIT2LTCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackST1T( int i ) const
// {
//   if( i>=0 && i<TrackST1TCol.size() )
//     return TrackST1TCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackST2T( int i ) const
// {
//   if( i>=0 && i<TrackST2TCol.size() )
//     return TrackST2TCol[i];
//   else
//     return 0;
// }

// inline TrLocalTrack * TrAnalyzer::GetTrackScatInT( int i ) const
// {
//   if( i>=0 && i<TrackScatInTCol.size() )
//     return TrackScatInTCol[i];
//   else
//     return 0;
// }

// inline BeamTrack * TrAnalyzer::GetBeamTrack( int i ) const
// {
//   if( i>=0 && i<BeamTrackCol.size() )
//     return BeamTrackCol[i];
//   else
//     return 0;
// }

// inline PreInTrack * TrAnalyzer::GetPreInTrack( int i ) const
// {
//   if( i>=0 && i<PreInTrackCol.size() )
//     return PreInTrackCol[i];
//   else
//     return 0;
// }

// inline PreOutTrack * TrAnalyzer::GetPreOutTrack( int i ) const
// {
//   if( i>=0 && i<PreOutTrackCol.size() )
//     return PreOutTrackCol[i];
//   else
//     return 0;
// }

// inline PreOut2Track * TrAnalyzer::GetPreOut2Track( int i ) const
// {
//   if( i>=0 && i<PreOut2TrackCol.size() )
//     return PreOut2TrackCol[i];
//   else
//     return 0;
// }

// inline Scat1Track * TrAnalyzer::GetScat1Track( int i ) const
// {
//   if( i>=0 && i<Scat1TrackCol.size() )
//     return Scat1TrackCol[i];
//   else
//     return 0;
// }

// inline Scat2Track * TrAnalyzer::GetScat2Track( int i ) const
// {
//   if( i>=0 && i<Scat2TrackCol.size() )
//     return Scat2TrackCol[i];
//   else
//     return 0;
// }

// inline Scat2ATrack * TrAnalyzer::GetScat2ATrack( int i ) const
// {
//   if( i>=0 && i<Scat2ATrackCol.size() )
//     return Scat2ATrackCol[i];
//   else
//     return 0;
// }

// inline Scat2BTrack * TrAnalyzer::GetScat2BTrack( int i ) const
// {
//   if( i>=0 && i<Scat2BTrackCol.size() )
//     return Scat2BTrackCol[i];
//   else
//     return 0;
// }

// inline Scat3Track * TrAnalyzer::GetScat3Track( int i ) const
// {
//   if( i>=0 && i<Scat3TrackCol.size() )
//     return Scat3TrackCol[i];
//   else
//     return 0;
// }

#endif 
