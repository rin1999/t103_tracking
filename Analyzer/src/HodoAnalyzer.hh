/*
  HodoAnalyzer.hh

  2024/04 K. Shirotori
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <vector>

#include "DetectorInfo.hh"
#include "RawData.hh"

class RawData;
class Hodo2Hit;
class Hodo1Hit;
class HodoCluster;

typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <Hodo1Hit*> Hodo1HitContainer;

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
   Hodo2HitContainer UTOFCont;
   Hodo2HitContainer DTOFCont;
   Hodo2HitContainer LTOFCont;
   Hodo2HitContainer T0Cont;
   Hodo2HitContainer T0rCont;
   Hodo1HitContainer BrefCont;
   Hodo2HitContainer T1Cont;
   Hodo1HitContainer BHTCont;
   
   HodoClusterContainer T0ClCont;
   
public:
   bool DecodeRawHits(RawData* rawData);
   
   bool DecodeUTOFHits(RawData* rawData);  
   bool DecodeDTOFHits(RawData* rawData);  
   bool DecodeLTOFHits(RawData* rawData);  
   bool DecodeT0Hits(RawData* rawData);  
   bool DecodeT0rHits(RawData* rawData);  
   bool DecodeBrefHits(RawData* rawData);  
   bool DecodeT1Hits(RawData* rawData);  
   bool DecodeBHTHits(RawData* rawData);  
   
   int GetNHitsUTOF( void ) const { return UTOFCont.size(); };
   int GetNHitsDTOF( void ) const { return DTOFCont.size(); };
   int GetNHitsLTOF( void ) const { return LTOFCont.size(); };
   int GetNHitsT0( void ) const { return T0Cont.size(); };
   int GetNHitsT0r( void ) const { return T0rCont.size(); };
   int GetNHitsBref( void ) const { return BrefCont.size(); };
   int GetNHitsT1( void ) const { return T1Cont.size(); };
   int GetNHitsBHT( void ) const { return BHTCont.size(); };
   
   inline Hodo2Hit * GetHitUTOF( int i ) const;
   inline Hodo2Hit * GetHitDTOF( int i ) const;
   inline Hodo2Hit * GetHitLTOF( int i ) const;
   inline Hodo2Hit * GetHitT0( int i ) const;
   inline Hodo2Hit * GetHitT0r( int i ) const;
   inline Hodo1Hit * GetHitBref( int i ) const;
   inline Hodo2Hit * GetHitT1( int i ) const;
   inline Hodo1Hit * GetHitBHT( int i ) const;
   
   int GetNClustersT0( void ) const { return T0ClCont.size(); }
   
   inline HodoCluster * GetClusterT0( int i )  const;
   
   bool ReCalcUTOFHits( bool applyRecursively=false );
   bool ReCalcDTOFHits( bool applyRecursively=false );
   bool ReCalcLTOFHits( bool applyRecursively=false );
   bool ReCalcT0Hits( bool applyRecursively=false );
   bool ReCalcT0rHits( bool applyRecursively=false );
   bool ReCalcBrefHits( bool applyRecursively=false );
   bool ReCalcT1Hits( bool applyRecursively=false );
   bool ReCalcBHTHits( bool applyRecursively=false );
   
   bool ReCalcT0Clusters( bool applyRecursively=false );
   
   bool ReCalcAll( void );
   
private:
   void clearUTOFHits();
   void clearDTOFHits();
   void clearLTOFHits();
   void clearT0Hits();
   void clearT0rHits( void );
   void clearBrefHits( void );
   void clearT1Hits( void );
   void clearBHTHits( void );

   static int MakeUpClusters( const Hodo2HitContainer & HitCont,
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

inline Hodo2Hit * HodoAnalyzer::GetHitUTOF( int i ) const
{
   if( i>=0 && i<UTOFCont.size() )
      return UTOFCont[i];
   else
      return 0;
}


inline Hodo2Hit * HodoAnalyzer::GetHitDTOF( int i ) const
{
   if( i>=0 && i<DTOFCont.size() )
      return DTOFCont[i];
   else
      return 0;
}


inline Hodo2Hit * HodoAnalyzer::GetHitLTOF( int i ) const
{
   if( i>=0 && i<LTOFCont.size() )
      return LTOFCont[i];
   else
      return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitT0( int i ) const
{
   if( i>=0 && i<T0Cont.size() )
      return T0Cont[i];
   else
      return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitT0r( int i ) const
{
   if( i>=0 && i<T0rCont.size() )
      return T0rCont[i];
   else
      return 0;
}

inline Hodo1Hit * HodoAnalyzer::GetHitBref( int i ) const
{
   if( i>=0 && i<BrefCont.size() )
      return BrefCont[i];
   else
      return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitT1( int i ) const
{
   if( i>=0 && i<T1Cont.size() )
      return T1Cont[i];
   else
      return 0;
}

inline Hodo1Hit * HodoAnalyzer::GetHitBHT( int i ) const
{
   if( i>=0 && i<BHTCont.size() )
      return BHTCont[i];
   else
      return 0;
}

#endif
