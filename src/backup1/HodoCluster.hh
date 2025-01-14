/*
  HodoCluster.hh

  2012/5  K.Shirotori
*/

#ifndef HODOCLUSTER_H
#define HODOCLUSTER_H

#include <cstddef>

class HodoHit;
class HodoAnalyzer;

class HodoCluster
{
public:
  HodoCluster( HodoHit *hitA, HodoHit *hitB=0, HodoHit *hitC=0 );
  virtual ~HodoCluster() {};

private:
  HodoCluster( const HodoCluster & );
  HodoCluster & operator = ( const HodoCluster & );

private:
  HodoHit *hitA_, *hitB_, *hitC_;
  int csize_;
  double MeanTime_, dE_;
  double MeanSeg_;
  double Tdiff_;
  bool gfastatus_;

public:
  int ClusterSize( void ) const { return csize_; };
  double CMeanTime( void ) const { return MeanTime_; };
  double DeltaE( void ) const { return dE_; };
  double MeanSeg( void ) const { return MeanSeg_; };
  double TimeDif( void ) const { return Tdiff_; };
  
  HodoHit * GetHit( int i ) const;
  bool GoodForAnalysis( void ) const { return gfastatus_; };
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } ;
  
  bool ReCalc( bool applyRecusively=false );
  
private:
  void calculate( void );
};
#endif
