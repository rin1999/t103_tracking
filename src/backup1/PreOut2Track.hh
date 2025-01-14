/*
  PreOut2Track.hh

  2016/2  K.Shirotori
*/

#ifndef PreOut2Track_h
#define PreOut2Track_h 1

#include "RungeKuttaUtilitiesPreOut2.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class PreOut2Track
{
public:
  PreOut2Track( double IniP,
		TrLocalTrack *Out,
		int type );
  ~PreOut2Track();

private:
  PreOut2Track( const PreOut2Track & );
  PreOut2Track & operator= ( const PreOut2Track & );

public:
  TrLocalTrack *GetLocalTrackOut( void ) { return Out_;}

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
  TrLocalTrack *Out_;
  double iniP_;
  int type_;

  std::vector <TrackHit *> hitArray_;
  RKPreOut2HitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  ThreeVector TofPos_, TofMom_;
  double pLenTtoT_, pLenTot_;
  RKPreOut2CordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKPreOut2HitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKPreOut2HitPointContainer &hpCont,
			    RKPreOut2CordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;
  
  void saveCalcPosition( const RKPreOut2HitPointContainer &hpCont );
  void PrintCalcHits( const RKPreOut2HitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKPreOut2CordParameter &cp );  
  bool doFitAgain( RKPreOut2CordParameter iniCord );
};

struct PreOut2TrackComp
  : public std::binary_function <PreOut2Track *, PreOut2Track *, bool>
{
  bool operator()( const PreOut2Track * const p1,
		   const PreOut2Track * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
