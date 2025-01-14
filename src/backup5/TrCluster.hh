/*
  TrCluster.hh

  2018/2  K.Shirotori
*/

#ifndef TrCluster_h
#define TrCluster_h

#include <cstddef>

class TrHit;
class TrAnalyzer;

class TrCluster
{
public:
  TrCluster( TrHit *hitA, TrHit *hitB=0, TrHit *hitC=0 );
  virtual ~TrCluster() {};
  
private:
  TrCluster( const TrCluster & );
  TrCluster & operator = ( const TrCluster & );
  
private:
  TrHit *hitA_, *hitB_, *hitC_;
  int csize_;
  double MeanTime_;
  double MeanFiber_;
  bool gfastatus_;

public:
  int ClusterSize( void ) const { return csize_; };
  double CMeanTime( void ) const { return MeanTime_; };
  double CMeanFiber( void ) const { return MeanFiber_; };
  
  TrHit * GetHit( int i ) const;
  bool GoodForAnalysis( void ) const { return gfastatus_; };
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } ;
  
  bool ReCalc( bool applyRecusively=false );
  
private:
  void calculate( void );
};

#endif
