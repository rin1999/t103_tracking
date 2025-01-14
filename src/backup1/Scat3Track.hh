/*
  Scat3Track.hh

  2016/2  K.Shirotori
*/

#ifndef Scat3Track_h
#define Scat3Track_h 1

#include "RungeKuttaUtilitiesScat3.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class Scat3Track
{
public:
  Scat3Track( double IniP,
  	      TrLocalTrack *In1, 
  	      TrLocalTrack *In2 );
  ~Scat3Track();

private:
  Scat3Track( const Scat3Track & );
  Scat3Track & operator= ( const Scat3Track & );

public:
  TrLocalTrack *GetLocalTrackIn1( void ) { return In1_;}
  TrLocalTrack *GetLocalTrackIn2( void ) { return In2_;}

  bool doFit( void );
  bool Status( void ) const { return status_; }
  int Niteration( void ) const { return nIteration_; }
  void SetInitialMomentum( double Pini ) { iniP_=Pini; }

  const ThreeVector &
  PrimaryPosition( void ) const { return priPos_; }
  const ThreeVector &
  PrimaryMomentum( void ) const { return priMom_; }
  double PathLengthToTOF( void ) const { return pLenTtoT_; }
  double PathLengthTotal( void ) const { return pLenTot_; }
  double chisqr( void ) const { return chisqr_; }

  std::size_t GetNHits( void ) const { return hitArray_.size(); } 
  TrackHit * GetHit( std::size_t nth ) const;
  TrackHit * GetHitOfLayerNumber( int lnum ) const;
  bool GoodForAnalysis( void ) const { return gfastatus_; }
  bool GoodForAnalysis( bool status )
  { bool ret=gfastatus_; gfastatus_=status; return ret; } 

  bool GetTrajectoryLocalPosition( int layer, double & x, double & y ) const;

  bool ReCalc( bool applyRecursively=false );

private:
  TrLocalTrack *In1_, *In2_;
  double iniP_;

  std::vector <TrackHit *> hitArray_;
  RKScat3HitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  double pLenTtoT_, pLenTot_;
  RKScat3CordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKScat3HitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKScat3HitPointContainer &hpCont,
			    RKScat3CordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;

  void saveCalcPosition( const RKScat3HitPointContainer &hpCont );
  void PrintCalcHits( const RKScat3HitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKScat3CordParameter &cp );  
  bool doFitAgain( RKScat3CordParameter iniCord );
};

struct Scat3TrackComp
  : public std::binary_function <Scat3Track *, Scat3Track *, bool>
{
  bool operator()( const Scat3Track * const p1,
		   const Scat3Track * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
