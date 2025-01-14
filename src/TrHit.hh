/*
  TrHit.hh

  2024/11  K.Shirotori
*/

#ifndef TrHit_h 
#define TrHit_h

#include <vector>
#include <deque>
#include <string>
#include <iostream>

typedef std::vector <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class TrLTrackHit;

class TrHit
{
public:
  TrHit();
  TrHit( int layer, double fiber ); 

  ~TrHit();

private:
  TrHit( const TrHit & );
  TrHit & operator = ( const TrHit & );

private:
  int layer_;
  int fiber_;
  double mfiber_;
  IntVec tdc_;
  IntVec tot_;
  IntVec trailing_;
  DoubleVec time_, dl_;
  DoubleVec trailingTime_;
  DoubleVec mtime_;
  double fpos_, mfpos_;
  double angle_;
  BoolVec belongTrack_;
  BoolVec timeRange_;
  int csize_;
  
  mutable std::vector <TrLTrackHit *> Cont_;
  
public:
  bool CalcTrObservables( void );
  bool ReCalcTrObservables( void );
  
  void SetLayer( int layer ) { layer_=layer; }
  void SetFiber( int fiber ) { fiber_=fiber; }  
  void SetTdcVal( int tdc );
  void SetTdcTrailing(int tdc) { trailing_.push_back(tdc); }
  void SetTotVal( int tot );

  void SetTime( double time ) { time_.push_back(time); }
  void SetDriftLength( double dl ) { dl_.push_back(dl); }
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  void SetTrailingTime( double t ) { trailingTime_.push_back(t); }

  void SetMeanFiber( double fiber ) { mfiber_=fiber; }
  void SetMeanTime( double mtime ) { mtime_.push_back(mtime); }
  void SetClusterSize( int size ) { csize_=size; }  

  int GetLayer( void ) const { return layer_; }
  int GetFiber( void )  const { return fiber_;  }
  double GetMeanFiber( void ) const { return mfiber_; }
  
  int GetTdcSize( void ) const { return tdc_.size(); }
  int GetTimeSize( void ) const { return time_.size(); }
  int GetDriftLengthSize( void ) const { return dl_.size(); }
  int GetTdcVal( int nh=0 ) const { return tdc_[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return trailing_[nh]; }
  int GetTotVal( int nh=0 ) const { return tot_[nh]; }

  double GetTime( int nh=0 ) const { return time_[nh]; }
  double GetDriftLength( int nh=0 ) const { return dl_[nh]; }
  double GetTrailingTime( int nh=0 ) const { return trailingTime_[nh]; }

  int GetMeanTimeSize( void ) const { return mtime_.size(); }
  double GetMeanTime( int nh=0 ) const { return mtime_[nh]; }
  int GetClusterSize( void ) const { return csize_; }

  double GetTiltAngle( void ) const { return angle_; }
  double GetPosition( void ) const { return fpos_; }
  double GetMPosition( void ) const { return mfpos_; }
  
  void setFlags( int nh=0 ) { belongTrack_[nh]=true; }
  void clearFlags( int nh=0 ) { belongTrack_[nh]=false; }
  bool showFlags( int nh=0 ) const { return belongTrack_[nh]; }
  bool rangecheck( int nh=0 ) const { return timeRange_[nh]; }
  void setRangeCheckStatus( bool status, 
			    int nh ) 
  { timeRange_[nh] = (status); }
  
  void RegisterHits( TrLTrackHit *hit ) const
  { Cont_.push_back(hit); }
  
  bool ReCalcTr( bool applyRecursively=false ) 
  { return CalcTrObservables(); }
  
private:
  void clearRegisteredHits( void );
  
};

#endif
