/*
  HodoAnalyzer.hh

  2012/5  K.Shirotori
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "s_BeamRawHit.hh"
#include "s_ScatRawHit.hh"

class RawData;
class HodoHit;
class HodoCluster;

typedef std::vector <HodoHit*> HodoHitContainer;
typedef std::vector <HodoCluster*> HodoClusterContainer;

class HodoAnalyzer
{
public:
  HodoAnalyzer();
  ~HodoAnalyzer();
private:
  HodoAnalyzer(const HodoAnalyzer &);
  HodoAnalyzer & operator =(const HodoAnalyzer &);

private:
  HodoHitContainer T0Cont;
  HodoHitContainer TOFCont;
  HodoHitContainer ITOFCont;
  HodoHitContainer PADCont;
  HodoHitContainer RICHCont;
  HodoHitContainer PID1Cont;
  HodoHitContainer PID2Cont;
  HodoHitContainer MFCont;
  HodoHitContainer VDCont;

  HodoClusterContainer T0ClCont;
  HodoClusterContainer TOFClCont;
  HodoClusterContainer ITOFClCont;
  HodoClusterContainer PADClCont;

public:
  bool DecodeRawHits(RawData* rawData);

  bool DecodeT0Hits(RawData* rawData);  
  bool DecodeTOFHits(RawData* rawData);  
  bool DecodeITOFHits(RawData* rawData);
  bool DecodePADHits(RawData* rawData);  
  bool DecodeRICHHits(RawData* rawData);  
  bool DecodePID1Hits(RawData* rawData);  
  bool DecodePID2Hits(RawData* rawData);  
  bool DecodeMFHits(RawData* rawData);  
  bool DecodeVDHits(RawData* rawData);  

  int GetNHitsT0( void ) const { return T0Cont.size(); };
  int GetNHitsTOF( void ) const { return TOFCont.size(); };
  int GetNHitsITOF( void ) const { return ITOFCont.size(); };
  int GetNHitsPAD( void ) const { return PADCont.size(); };
  int GetNHitsRICH( void ) const { return RICHCont.size(); };
  int GetNHitsPID1( void ) const { return PID1Cont.size(); };
  int GetNHitsPID2( void ) const { return PID2Cont.size(); };
  int GetNHitsMF( void ) const { return MFCont.size(); };
  int GetNHitsVD( void ) const { return VDCont.size(); };

  inline HodoHit * GetHitT0( int i ) const;
  inline HodoHit * GetHitTOF( int i ) const;
  inline HodoHit * GetHitITOF( int i ) const;
  inline HodoHit * GetHitPAD( int i ) const;
  inline HodoHit * GetHitRICH( int i ) const;
  inline HodoHit * GetHitPID1( int i ) const;
  inline HodoHit * GetHitPID2( int i ) const;
  inline HodoHit * GetHitMF( int i ) const;
  inline HodoHit * GetHitVD( int i ) const;

  int GetNClustersT1( void ) const { return T0ClCont.size(); };
  int GetNClustersTOF( void ) const { return TOFClCont.size(); };
  int GetNClustersITOF( void ) const { return ITOFClCont.size(); };
  int GetNClustersPAD( void ) const { return PADClCont.size(); };

  inline HodoCluster * GetClusterT0( int i ) const;
  inline HodoCluster * GetClusterTOF( int i ) const;
  inline HodoCluster * GetClusterITOF( int i ) const;
  inline HodoCluster * GetClusterPAD( int i ) const;

  bool ReCalcT0Hits( bool applyRecursively=false );
  bool ReCalcTOFHits( bool applyRecursively=false );
  bool ReCalcITOFHits( bool applyRecursively=false );
  bool ReCalcPADHits( bool applyRecursively=false );
  bool ReCalcRICHHits( bool applyRecursively=false );
  bool ReCalcPID1Hits( bool applyRecursively=false );
  bool ReCalcPID2Hits( bool applyRecursively=false );
  bool ReCalcMFHits( bool applyRecursively=false );
  bool ReCalcVDHits( bool applyRecursively=false );
  
  bool ReCalcT0Clusters( bool applyRecursively=false );
  bool ReCalcTOFClusters( bool applyRecursively=false );
  bool ReCalcITOFClusters( bool applyRecursively=false );
  bool ReCalcPADClusters( bool applyRecursively=false );

  bool ReCalcAll( void );

private:
  void clearT0Hits();
  void clearTOFHits();
  void clearITOFHits();
  void clearPADHits();
  void clearRICHHits();
  void clearPID1Hits();
  void clearPID2Hits();
  void clearMFHits();
  void clearVDHits();

  static int MakeUpClusters( const HodoHitContainer & HitCont,
			     HodoClusterContainer & ClusterCont,
			     double maxTimeDif );
};

inline HodoCluster * HodoAnalyzer::GetClusterT0( int i ) const
{
  if( i>=0 && i<T0ClCont.size() )
    return T0ClCont[i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterTOF( int i ) const
{
  if( i>=0 && i<TOFClCont.size() )
    return TOFClCont[i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterITOF( int i ) const
{
  if( i>=0 && i<ITOFClCont.size() )
    return ITOFClCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitT0( int i ) const
{
  if( i>=0 && i<T0Cont.size() )
    return T0Cont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitTOF( int i ) const
{
  if( i>=0 && i<TOFCont.size() )
    return TOFCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitITOF( int i ) const
{
  if( i>=0 && i<ITOFCont.size() )
    return ITOFCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitPAD( int i ) const
{
  if( i>=0 && i<PADCont.size() )
    return PADCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitRICH( int i ) const
{
  if( i>=0 && i<RICHCont.size() )
    return RICHCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitPID1( int i ) const
{
  if( i>=0 && i<PID1Cont.size() )
    return PID1Cont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitPID2( int i ) const
{
  if( i>=0 && i<PID2Cont.size() )
    return PID2Cont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitMF( int i ) const
{
  if( i>=0 && i<MFCont.size() )
    return MFCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitVD( int i ) const
{
  if( i>=0 && i<VDCont.size() )
    return VDCont[i];
  else
    return 0;
}

#endif
