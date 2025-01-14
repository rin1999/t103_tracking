/*
  DCHit.hh

  2018/12  K.Shirotori
*/

#ifndef DCHit_h 
#define DCHit_h

#include <vector>
#include <deque>
#include <string>
#include <iostream>

typedef std::vector <bool> BoolVec;
typedef std::vector <int> IntVec;
typedef std::vector <double> DoubleVec;

class DCLTrackHit;

class DCHit
{
public:
  DCHit();
  DCHit( int layer, double wire ); 

  ~DCHit();

private:
  DCHit( const DCHit & );
  DCHit & operator = ( const DCHit & );

private:
  int layer_;
  double wire_;
  IntVec tdc_;
  IntVec trailing_;
  DoubleVec dt_, dl_;
  DoubleVec trailingTime_;
  
  double wpos_;
  double angle_;
  BoolVec belongTrack_;
  BoolVec dlRange_;
  
  mutable std::vector <DCLTrackHit *> Cont_;
  
public:
  bool CalcDCObservables( void );
  
  void SetLayer( int layer ) { layer_=layer; }
  void SetWire( double wire ) { wire_=wire; }  
  void SetTdcVal( int tdc );
  void SetTdcTrailing(int tdc) {trailing_.push_back(tdc); }
  void SetDriftTime( double dt ) { dt_.push_back(dt); }
  void SetDriftLength( double dl ) { dl_.push_back(dl); }
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  void SetTrailingTime( double t ) { trailingTime_.push_back(t); }
  
  int GetLayer( void ) const { return layer_; }
  double GetWire( void )  const { return int(wire_);  }
  
  int GetTdcSize( void ) const { return tdc_.size(); }
  int GetDriftTimeSize( void ) const { return dt_.size(); }
  int GetDriftLengthSize( void ) const { return dl_.size(); }
  int GetTdcVal( int nh=0 ) const { return tdc_[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return trailing_[nh]; }
  double GetDriftTime( int nh=0 ) const { return dt_[nh]; }
  double GetDriftLength( int nh=0 ) const { return dl_[nh]; }
  double GetTrailingTime( int nh=0 ) const { return trailingTime_[nh]; }
  
  double GetTiltAngle( void ) const { return angle_; }
  double GetWirePosition( void ) const { return wpos_; }
  
  void setFlags( int nh=0 ) { belongTrack_[nh]=true; }
  void clearFlags( int nh=0 ) { belongTrack_[nh]=false; }
  bool showFlags( int nh=0 ) const { return belongTrack_[nh]; }
  bool rangecheck( int nh=0 ) const { return dlRange_[nh]; }
  void setRangeCheckStatus( bool status, 
			    int nh ) 
  { dlRange_[nh] = (status); }
  
  void RegisterHits( DCLTrackHit *hit ) const
  { Cont_.push_back(hit); }
  
  bool ReCalcDC( bool applyRecursively=false ) 
  { return CalcDCObservables(); }
  
private:
  void clearRegisteredHits( void );
  
};

#endif
