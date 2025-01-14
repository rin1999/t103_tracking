/*
  PreOutTrack.hh

  2016/2  K.Shirotori
*/

#ifndef PreOutTrack_h
#define PreOutTrack_h 1

#include "RungeKuttaUtilitiesPreOut.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class PreOutTrack
{
public:
  PreOutTrack( double IniP,
	       TrLocalTrack *Out1, 
	       TrLocalTrack *Out2 );
  ~PreOutTrack();

private:
  PreOutTrack( const PreOutTrack & );
  PreOutTrack & operator= ( const PreOutTrack & );

public:
  TrLocalTrack *GetLocalTrackOut1( void ) { return Out1_;}
  TrLocalTrack *GetLocalTrackOut2( void ) { return Out2_;}

  bool doFit( void );
  bool Status( void ) const { return status_; }
  int Niteration( void ) const { return nIteration_; }
  void SetInitialMomentum( double Pini ) { iniP_=Pini; }

  const ThreeVector &
  PrimaryPosition( void ) const { return priPos_; }
  const ThreeVector &
  PrimaryMomentum( void ) const { return priMom_; }
  const ThreeVector &
  TofPosition( void ) const { return TofPos_; }
  const ThreeVector &
  TofMomentum( void ) const { return TofMom_; }
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
  TrLocalTrack *Out1_, *Out2_;
  double iniP_;

  std::vector <TrackHit *> hitArray_;
  RKPreOutHitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;  
  ThreeVector TofPos_, TofMom_;
  double pLenTtoT_, pLenTot_;
  RKPreOutCordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKPreOutHitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKPreOutHitPointContainer &hpCont,
			    RKPreOutCordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;
  
  void saveCalcPosition( const RKPreOutHitPointContainer &hpCont );
  void PrintCalcHits( const RKPreOutHitPointContainer &hpCont,
		      std::ostream &ost ) const;
  
  bool saveTrackParameters( const RKPreOutCordParameter &cp );  
  bool doFitAgain( RKPreOutCordParameter iniCord );
};

struct PreOutTrackComp
  : public std::binary_function <PreOutTrack *, PreOutTrack *, bool>
{
  bool operator()( const PreOutTrack * const p1,
		   const PreOutTrack * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
