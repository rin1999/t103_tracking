/*
  Scat2Track.hh

  2016/2  K.Shirotori
*/

#ifndef Scat2Track_h
#define Scat2Track_h 1

#include "RungeKuttaUtilitiesScat2.hh"
#include "ThreeVector.hh"

#include <vector>
#include <iosfwd>
#include <functional>

class PreInTrack;
class PreOut2Track;
class TrLocalTrack;
class TrackHit;
class TrAnalyzer;

class Scat2Track
{
public:
  Scat2Track( double IniP,
	      PreInTrack *In, 
	      PreOut2Track *Out,
	      int type );
  ~Scat2Track();

private:
  Scat2Track( const Scat2Track & );
  Scat2Track & operator= ( const Scat2Track & );

public:
  PreInTrack *GetPreInTrack( void ) { return In_;}
  PreOut2Track *GetPreOut2TrackOut( void ) { return Out_;}

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
  PreInTrack *In_;
  PreOut2Track *Out_;
  double iniP_;
  int type_;

  std::vector <TrackHit *> hitArray_;
  RKScat2HitPointContainer hPntCont_;

  bool status_;
  int nIteration_;
  double chisqr_;

  ThreeVector priPos_, priMom_;
  double pLenTtoT_, pLenTot_;
  RKScat2CordParameter CPval_;

  bool gfastatus_;

private:
  void fillHitArray( void );
  void clearHitArray( void );
  double calcChiSqr( const RKScat2HitPointContainer &hpCont ) const;
  bool guessNextParameters( const RKScat2HitPointContainer &hpCont,
			    RKScat2CordParameter &Cord,
			    double &estDeltaChisqr,
			    double &lambdaCri, double dmp=0.0 ) const;

  void saveCalcPosition( const RKScat2HitPointContainer &hpCont );
  void PrintCalcHits( const RKScat2HitPointContainer &hpCont,
		      std::ostream &ost ) const;

  bool saveTrackParameters( const RKScat2CordParameter &cp );  
  bool doFitAgain( RKScat2CordParameter iniCord );
};

struct Scat2TrackComp
  : public std::binary_function <Scat2Track *, Scat2Track *, bool>
{
  bool operator()( const Scat2Track * const p1,
		   const Scat2Track * const p2 ) const
  {
    int n1=p1->GetNHits(), n2=p2->GetNHits();
    if( n1>n2+4 ) return true;
    else if( n2>n1+4 ) return false;
    else
      return (p1->chisqr())<(p2->chisqr());
  }
};

#endif
