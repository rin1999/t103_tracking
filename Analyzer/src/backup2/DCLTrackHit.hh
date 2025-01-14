/*
  DCLTrackHit.hh

  2018/12  K.Shirotori
*/

#ifndef DCLTrackHit_h
#define DCLTrackHit_h 1

#include "DCHit.hh"

class DCAnalyzer;

class DCLTrackHit
{
public:
  DCLTrackHit( DCHit *hit, double pos, int nh ) 
    : Hit_(hit), xl_(pos), nth_Hit_(nh)
  { 
    Hit_->RegisterHits(this);
 }
  DCLTrackHit( const DCLTrackHit & right )
    :Hit_(right.Hit_), nth_Hit_(right.nth_Hit_), xl_(right.xl_), xcl_(right.xcl_), 
      xcal_(right.xcal_), ycal_(right.ycal_)
  { 
    Hit_->RegisterHits(this); 
}
private:
  ~DCLTrackHit() {
}

private:
  DCHit *Hit_;
  int    nth_Hit_;
  double xl_, xcl_, xcal_, ycal_;

public:
  void SetLocalHitPos( double xl ) { xl_=xl; }
  void SetCalPosition( double x, double y ) { xcal_=x; ycal_=y; }

  int GetLayer( void ) const { return Hit_->GetLayer(); }
  double GetWire( void ) const { return Hit_->GetWire(); }
  int GetTdcVal( void ) const { return Hit_->GetTdcVal(nth_Hit_); }
  int GetTdcSize( void ) const { return Hit_->GetTdcSize(); }

  double GetDriftTime( void ) const { return Hit_->GetDriftTime(nth_Hit_); }
  double GetDriftLength( void ) const { return Hit_->GetDriftLength(nth_Hit_); }
  double GetTiltAngle( void ) const { return Hit_->GetTiltAngle(); }
  double GetWirePosition( void ) const { return Hit_->GetWirePosition(); }
  double GetLocalHitPos( void ) const { return xl_; }
  double GetLocalCalPos( void ) const ;

  double GetXcal( void ) const { return xcal_; }
  double GetYcal( void ) const { return ycal_; }
  double GetResidual( void ) const { return xl_-GetLocalCalPos(); }

  void setFlags( void ) { Hit_->setFlags(nth_Hit_); }
  void clearFlags( void ) { Hit_->clearFlags(nth_Hit_); }
  bool showFlags( void ) const { return Hit_->showFlags(nth_Hit_); }

  bool ReCalc( bool applyRecursively=false ); 

  friend class DCHit;

};

#endif
