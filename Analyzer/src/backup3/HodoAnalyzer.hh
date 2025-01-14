/*
  HodoAnalyzer.hh
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class Hodo2Hit;
class Hodo1Hit;
class BH2Hit;
class HodoCluster;
class BH2Cluster;

typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <Hodo1Hit*> Hodo1HitContainer;
typedef std::vector <BH2Hit*>   BH2HitContainer;

typedef std::vector <HodoCluster*> HodoClusterContainer;
typedef std::vector <BH2Cluster*>  BH2ClusterContainer;

class HodoAnalyzer
{
public:
  HodoAnalyzer();
  ~HodoAnalyzer();
private:
  HodoAnalyzer(const HodoAnalyzer &);
  HodoAnalyzer & operator =(const HodoAnalyzer &);

private:
  Hodo1HitContainer GCCont;
  Hodo2HitContainer BH1Cont;
  BH2HitContainer   BH2Cont;
  Hodo1HitContainer BACCont;
  //  Hodo1HitContainer TgtCont;
  Hodo2HitContainer TOFCont;
  Hodo2HitContainer LCCont;
  Hodo1HitContainer ACCont[NumOfLayersAc+1];  

  HodoClusterContainer BH1ClCont;
  BH2ClusterContainer  BH2ClCont;
  HodoClusterContainer TOFClCont;
  HodoClusterContainer LCClCont;

public:
  bool DecodeRawHits(RawData* rawData);

  bool DecodeGCHits(RawData* rawData);  
  bool DecodeBH1Hits(RawData* rawData);
  bool DecodeBH2Hits(RawData* rawData);
  bool DecodeBACHits(RawData* rawData);
  //  bool DecodeTGTHits(RawData* rawData);
  bool DecodeTOFHits(RawData* rawData);
  bool DecodeLCHits(RawData* rawData);
  bool DecodeACHits(RawData* rawData);

  int GetNHitsGC( void ) const { return GCCont.size(); };
  int GetNHitsBH1( void ) const { return BH1Cont.size(); };
  int GetNHitsBH2( void ) const { return BH2Cont.size(); };
  int GetNHitsBAC( void ) const { return BACCont.size(); };
  int GetNHitsTOF( void ) const { return TOFCont.size(); };
  int GetNHitsLC( void ) const { return LCCont.size(); };
  inline int GetNHitsAC( int layer ) const;

  inline Hodo1Hit * GetHitGC( int i ) const;
  inline Hodo2Hit * GetHitBH1( int i ) const;
  inline BH2Hit   * GetHitBH2( int i ) const;
  inline Hodo1Hit * GetHitBAC( int i ) const;
  inline Hodo2Hit * GetHitTOF( int i ) const;
  inline Hodo2Hit * GetHitLC( int i ) const;
  inline Hodo1Hit * GetHitAC( int layer, int i ) const;

  int GetNClustersBH1( void ) const { return BH1ClCont.size(); };
  int GetNClustersBH2( void ) const { return BH2ClCont.size(); };
  int GetNClustersTOF( void ) const { return TOFClCont.size(); }
  int GetNClustersLC( void )  const { return LCClCont.size(); }

  inline HodoCluster * GetClusterBH1( int i ) const;
  inline BH2Cluster  * GetClusterBH2( int i ) const;
  inline HodoCluster * GetClusterTOF( int i ) const;
  inline HodoCluster * GetClusterLC( int i )  const;

  bool ReCalcGCHits( bool applyRecursively=false );
  bool ReCalcBH1Hits( bool applyRecursively=false );
  bool ReCalcBH2Hits( bool applyRecursively=false );
  bool ReCalcBACHits( bool applyRecursively=false );
  bool ReCalcTOFHits( bool applyRecursively=false );
  bool ReCalcLCHits( bool applyRecursively=false );
  bool ReCalcACHits( bool applyRecursively=false );
  
  bool ReCalcBH1Clusters( bool applyRecursively=false );
  bool ReCalcBH2Clusters( bool applyRecursively=false );
  bool ReCalcTOFClusters( bool applyRecursively=false );
  bool ReCalcLCClusters( bool applyRecursively=false );

  bool ReCalcAll( void );

private:
  void clearGCHits();
  void clearBH1Hits();
  void clearBH2Hits();
  void clearBACHits();
  void clearTOFHits( void );
  void clearLCHits( void );
  void clearACHits( void );

  static int MakeUpClusters( const Hodo2HitContainer & HitCont,
			     HodoClusterContainer & ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const BH2HitContainer & HitCont,
			     BH2ClusterContainer & ClusterCont,
			     double maxTimeDif );

};

inline HodoCluster * HodoAnalyzer::GetClusterBH1( int i ) const
{
  if( i>=0 && i<BH1ClCont.size() )
    return BH1ClCont[i];
  else
    return 0;
}

inline BH2Cluster * HodoAnalyzer::GetClusterBH2( int i ) const
{
  if( i>=0 && i<BH2ClCont.size() )
    return BH2ClCont[i];
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

inline HodoCluster * HodoAnalyzer::GetClusterLC( int i ) const
{
  if( i>=0 && i<LCClCont.size() )
    return LCClCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitBH1( int i ) const
{
  if( i>=0 && i<BH1Cont.size() )
    return BH1Cont[i];
  else
    return 0;
}

inline BH2Hit * HodoAnalyzer::GetHitBH2( int i ) const
{
  if( i>=0 && i<BH2Cont.size() )
    return BH2Cont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitTOF( int i ) const
{
  if( i>=0 && i<TOFCont.size() )
    return TOFCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitLC( int i ) const
{
  if( i>=0 && i<LCCont.size() )
    return LCCont[i];
  else
    return 0;
}

inline Hodo1Hit * HodoAnalyzer::GetHitAC( int layer, int i ) const
{
  if( layer==1 || layer==2 )
    {
      if( i>=0 && i<ACCont[layer].size() )
	return ACCont[layer][i];
      else
	return 0;
    }
  else 
    return 0;
}

inline int HodoAnalyzer::GetNHitsAC( int layer ) const
{
  if( layer==1 || layer==2 )
    return ACCont[layer].size();
  else
    return ACCont[0].size();
}


#endif
