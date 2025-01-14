/*
  PreInTrack.hh

  2016/2  K.Shirotori
*/

#ifndef PreInTrack_h
#define PreInTrack_h 1

#include "RungeKuttaUtilitiesPreIn.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class PreInTrack
{
public:
  PreInTrack( double IniP,
	      TrLocalTrack *In1, 
	      TrLocalTrack *In2 );
  ~PreInTrack();

private:
  PreInTrack( const PreInTrack & );
  PreInTrack & operator= ( const PreInTrack & );

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
  RKPreInHitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  double pLenTtoT_, pLenTot_;
  RKPreInCordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKPreInHitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKPreInHitPointContainer &hpCont,
			    RKPreInCordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;

  void saveCalcPosition( const RKPreInHitPointContainer &hpCont );
  void PrintCalcHits( const RKPreInHitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKPreInCordParameter &cp );  
  bool doFitAgain( RKPreInCordParameter iniCord );
};

struct PreInTrackComp
  : public std::binary_function <PreInTrack *, PreInTrack *, bool>
{
  bool operator()( const PreInTrack * const p1,
		   const PreInTrack * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
