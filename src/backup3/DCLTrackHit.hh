/*
  DCLTrackHit.hh
*/

#ifndef DCLTrackHit_h
#define DCLTrackHit_h 1

#include "DCHit.hh"

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class DCAnalyzer;

class DCLTrackHit
{
public:
  DCLTrackHit( DCHit *hit, double pos, int nh ) 
    : Hit_(hit), xl_(pos), nth_Hit_(nh)
  { 
    Hit_->RegisterHits(this);
#ifdef MemoryCheck
    ++sm_counter;
#endif 
 }
  DCLTrackHit( const DCLTrackHit & right )
    :Hit_(right.Hit_), nth_Hit_(right.nth_Hit_), xl_(right.xl_), xcl_(right.xcl_), 
      xcal_(right.xcal_), ycal_(right.ycal_)
  { 
    Hit_->RegisterHits(this); 
#ifdef MemoryCheck
    ++sm_counter;
#endif 
}
private:
  ~DCLTrackHit() {
#ifdef MemoryCheck
    --sm_counter;
#endif    
}

private:
  DCHit *Hit_;
  int    nth_Hit_;
  double xl_, xcl_, xcal_, ycal_;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  void SetLocalHitPos( double xl ) { xl_=xl; }
  //Add Y.Y.
  void SetLocalCalPosVXU( double xcl ) { xcl_=xcl; }
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
  //Add Y.Y.
  double GetLocalCalPosVXU( void ) const { return xcl_; }

  double GetXcal( void ) const { return xcal_; }
  double GetYcal( void ) const { return ycal_; }
  double GetResidual( void ) const { return xl_-GetLocalCalPos(); }
  //Add Y.Y.
  double GetResidualVXU( void ) const { return xl_-xcl_; }

  void setFlags( void ) { Hit_->setFlags(nth_Hit_); }
  void clearFlags( void ) { Hit_->clearFlags(nth_Hit_); }
  bool showFlags( void ) const { return Hit_->showFlags(nth_Hit_); }

  bool ReCalc( bool applyRecursively=false ); 

  friend class DCHit;

  // for DS
// private:
//   unsigned int DSKey_;
//   struct DSFormat {
//     unsigned int header_;
//     unsigned int pkey_;
//     float xl_, xcal_, ycal_;
//   };

// public:
//   DCLTrackHit( unsigned int *bufp, DCAnalyzer *DCana )
//   { DSRestore( bufp, DCana ); }
//   std::size_t DSSize( void ) const
//   {
//     return sizeof(DSFormat)/sizeof(int)+
//       ( sizeof(DSFormat)%sizeof(int) ? 1 : 0 );
//   }
//   std::size_t DSSave( unsigned int *bufp ) const;
//   bool DSRestore( unsigned int *bufp, DCAnalyzer *DCana );
//   void DSSetSaveKey( unsigned int key ) { DSKey_ = key; }
//   unsigned int DSGetSaveKey( void ) const { return DSKey_; } 
};

#endif
