/*
  DCLocalTrack.hh
*/

#ifndef DCLocalTrack_h
#define DCLocalTrack_h 1

#include <vector>
#include <functional>
#include "ThreeVector.hh"

class DCLTrackHit;
class DCAnalyzer;

class DCLocalTrack
{
public:
  explicit DCLocalTrack();
  ~DCLocalTrack();
private:
  DCLocalTrack( const DCLocalTrack & );
  DCLocalTrack & operator = ( const DCLocalTrack & );

private:
  std::vector <DCLTrackHit *> hitArray;

public:
  void AddHit( DCLTrackHit *hitp ) { hitArray.push_back( hitp ); }
  bool DoFit( void );  
  bool DoFitX( void );
  std::size_t GetNHit( void ) const { return hitArray.size(); }
  DCLTrackHit * GetHit( std::size_t nth ) const;
  DCLTrackHit * GetHitOfLayerNumber( int lnum ) const;

  double GetX0( void ) const { return x0_; }
  double GetY0( void ) const { return y0_; }
  double GetU0( void ) const { return u0_; }
  double GetV0( void ) const { return v0_; }

  double GetChiSquare( void ) const { return chisqr_; }
  double GetX( double z ) const { return x0_+u0_*z; } 
  double GetY( double z ) const { return y0_+v0_*z; } 
  bool GetStatus( void ) const { return status_; } 
  bool GoodForTracking( void ) const { return gftstatus_; }
  bool GoodForTracking( bool status )
  { bool ret=gftstatus_; gftstatus_=status; return ret; } 
  bool ReCalc( bool ApplyRecursively=false );  
private:
  bool status_;
  double x0_, y0_, u0_, v0_;
  double a_,b_;
  double chisqr_;
  bool gftstatus_;
};

struct DCLTrackComp 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    if( (n1>n2+1) ){
      return true;
    }
    else if( (n2>n1+1)  ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};

struct DCLTrackComp1 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    if(n1>n2) return true;
    else if(n2>n1) return false;
    else
      return (p1->GetChiSquare())<=(p2->GetChiSquare());
  }
};

struct DCLTrackComp2 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    if(n1<n2) return true;
    else if(n2<n1) return false;
    else
      return (p1->GetChiSquare())<=(p2->GetChiSquare());
  }
};

struct DCLTrackComp3 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    double a1=fabs(1.-chi1),a2=fabs(1.-chi2);
    if(a1<a2) return true;
    else if(a2<a1) return false;
    else
      return (n1<=n2);
  }
};

struct DCLTrackComp4 
  : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
{
  bool operator()( const DCLocalTrack * const p1, 
		   const DCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    //if( (n1>n2+1) ){
    //    if( (n1>n2+1) && (fabs(chi1-chi2)<5.) ){
    if( (n1>n2+1) && (fabs(chi1-chi2)<2.) ){
      return true;
    }
    //else if( (n2>n1+1)  ){
    //    else if( (n2>n1+1) && (fabs(chi1-chi2)<5.) ){
    else if( (n2>n1+1) && (fabs(chi1-chi2)<2.) ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};



#endif
