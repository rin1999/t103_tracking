/*
  TrackHit.hh
*/


#ifndef TrackHit_h
#define TrackHit_h 1

#include "DCLTrackHit.hh" 
#include "ThreeVector.hh"

#include <cstddef>

class DCLocalTrack;

class TrackHit
{
public:
  explicit TrackHit( DCLTrackHit *hit );
  ~TrackHit();
private:
  TrackHit( const TrackHit & );
  TrackHit & operator= ( const TrackHit & );


public:
  DCLTrackHit *GetHit( void ) { return dchitp_; }
  void SetCalGPos( const ThreeVector &pos ) { calGPos_=pos; }
  void SetCalLPos( double pos ) { calLPos_=pos; }

  int GetLayer( void ) const { return dchitp_->GetLayer(); }
  double GetLocalHitPos( void ) const 
  { return dchitp_->GetLocalHitPos(); }
  const ThreeVector & GetCalGPos( void ) const 
  { return calGPos_; }
  double GetCalLPos( void ) const { return calLPos_; }
  double GetResidual( void ) const 
  { return dchitp_->GetLocalHitPos()-calLPos_; }

  bool ReCalc( bool applyRecursively=false )
  {
    if( applyRecursively )
      return dchitp_->ReCalc(applyRecursively);
    else
      return true;
  }

private:
  DCLTrackHit *dchitp_;
  ThreeVector calGPos_;
  double calLPos_;

  // for DS
// private:
//   unsigned int DSKey_;
//   struct DSFormat {
//     unsigned int header_;
//     unsigned int pkey_;
//     float gxcal_, gycal_, gzcal_, callpos_;
//   };    

// public:
//   explicit TrackHit( unsigned int *bufp, 
// 		     DCLocalTrack *trIn, DCLocalTrack *trOut )
//   { DSRestore( bufp, trIn, trOut ); }

//   std::size_t DSSize( void ) const
//   {
//     return sizeof(DSFormat)/sizeof(int)+
//       ( sizeof(DSFormat)%sizeof(int) ? 1 : 0 );
//   }
//   std::size_t DSSave( unsigned int *bufp ) const;
//   bool DSRestore( unsigned int *bufp, 
// 		  DCLocalTrack *trIn, DCLocalTrack *trOut );
//   void DSSetSaveKey( unsigned int key ) { DSKey_ = key; }
//   unsigned int DSGetSaveKey( void ) const { return DSKey_; } 
};

#endif

