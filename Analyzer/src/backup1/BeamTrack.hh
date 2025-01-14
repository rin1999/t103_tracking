/*
  BeamTrack.hh

  2016/2  K.Shirotori
*/

#ifndef BeamTrack_h
#define BeamTrack_h 1

#include "RungeKuttaUtilitiesBeam.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class BeamTrack
{
public:
  BeamTrack( TrLocalTrack *In );
  ~BeamTrack();

private:
  BeamTrack( const BeamTrack & );
  BeamTrack & operator= ( const BeamTrack & );

public:
  TrLocalTrack *GetLocalTrackIn( void ) { return In_;}
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
  TrLocalTrack *In_;
  double iniP_;

  std::vector <TrackHit *> hitArray_;
  RKBeamHitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  double pLenTtoT_, pLenTot_;
  RKBeamCordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKBeamHitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKBeamHitPointContainer &hpCont,
			    RKBeamCordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;

  void saveCalcPosition( const RKBeamHitPointContainer &hpCont );
  void PrintCalcHits( const RKBeamHitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKBeamCordParameter &cp );  
  bool doFitAgain( RKBeamCordParameter iniCord );
};

struct BeamTrackComp
  : public std::binary_function <BeamTrack *, BeamTrack *, bool>
{
  bool operator()( const BeamTrack * const p1,
		   const BeamTrack * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
