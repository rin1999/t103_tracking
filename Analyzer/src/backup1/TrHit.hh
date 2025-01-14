/*
  TrHit.hh

  2016/2  K.Shirotori
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
  TrHit( int layer, double wire ); 

  ~TrHit();

private:
  TrHit( const TrHit & );
  TrHit & operator = ( const TrHit & );

private:
  int layer_;
  double wire_;
  DoubleVec pos_;
  DoubleVec dt_, dl_;
//   IntVec tdc_;
//   IntVec trailing_;
//   DoubleVec dt_, dl_;
//   DoubleVec trailingTime_;

  double wpos_;
  double angle_;
  BoolVec belongTrack_;
  BoolVec dlRange_;

//   int clsize_;
//   bool mwpcflag_;
//   double mwire_;
//   double mwpos_;

  mutable std::vector <TrLTrackHit *> Cont_;

public:
  bool CalcObservables( void );

  void SetLayer( int layer ) { layer_=layer; }
  void SetWire( double wire ) { wire_=wire; }  
  void SetPos( double pos );

  void SetTdcVal( int tdc );
  //   void SetTdcTrailing(int tdc) {trailing_.push_back(tdc); }
  void SetDriftTime( double dt ) { dt_.push_back(0.); }
  void SetDriftLength( double dl ) { dl_.push_back(dl); }
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  //  void SetTrailingTime( double t ) { trailingTime_.push_back(t); }

//   void SetClusterSize( int size ) { clsize_=size; }
//   void SetMWPCFlag( bool flag ) { mwpcflag_=flag; }
//   void SetMeanWire( double mwire ) { mwire_=mwire; }
//   void SetMeanWirePosition( double mwpos ) { mwpos_=mwpos; }
//   void SetWirePosition( double wpos ) { wpos_ = wpos; }

  int GetLayer( void ) const { return layer_; }
  double GetWire( void )  const { return wire_;}

  int GetPosSize( void ) const { return pos_.size(); }
  int GetDriftTimeSize( void ) const { return pos_.size(); }
  int GetDriftLengthSize( void ) const { return pos_.size(); }
  double GetPos( int nh=0 ) const { return pos_[nh]; }
//   int GetTdcTrailing( int nh=0 ) const { return trailing_[nh]; }
  double GetDriftTime( int nh=0 ) const { return 0.; }
  double GetDriftLength( int nh=0 ) const { return dl_[nh]; }
//   double GetTrailingTime( int nh=0 ) const { return trailingTime_[nh]; }

  double GetTiltAngle( void ) const { return angle_; }
  double GetWirePosition( void ) const { return wpos_; }

//   int GetClusterSize( void ) const { return clsize_; }
//   double GetMeamWire( void ) const { return mwire_; }
//   double GetMeamWirePosition( void ) const { return mwpos_; }

  void setFlags( int nh=0 ) { belongTrack_[nh]=true; }
  void clearFlags( int nh=0 ) { belongTrack_[nh]=false; }
  bool showFlags( int nh=0 ) const { return belongTrack_[nh]; }
  bool rangecheck( int nh=0 ) const { return dlRange_[nh]; }
  void setRangeCheckStatus( bool status, 
 			    int nh ) 
  { dlRange_[nh] = (status); }

  void RegisterHits( TrLTrackHit *hit ) const
  { Cont_.push_back(hit); }

  bool ReCalc( bool applyRecursively=false ) 
  { return CalcObservables(); }
//   bool ReCalcMWPC( bool applyRecursively=false ) 
//   { return CalcMWPCObservables(); }

private:
  void clearRegisteredHits( void );

};


#endif
