/*
  RungeKuttaUtilitiesBeam.hh

  2004/6/28    T.Takahashi

*/

#ifndef RungeKuttaUtilitiesBeam_h

#define RungeKuttaUtilitiesBeam_h 1

#include "ThreeVector.hh"

#include <vector>
#include <utility>
#include <iosfwd>

class RKBeamFieldIntegral;
class RKBeamDeltaFieldIntegral;
class RKBeamTrajectoryPoint;
class RKBeamcalcHitPoint;
class RKBeamCordParameter;
class RKBeamHitPointContainer;

class RKBeamFieldIntegral
{
public:
  RKBeamFieldIntegral( double Ky, double Kz, 
		       double Ayu, double Ayv, double Azu, double Azv,
		       double Cyy=0.0, double Cyz=0.0, 
		       double Czy=0.0, double Czz=0.0 )
    : ky(Ky), kz(Kz), ayu(Ayu), ayv(Ayv), azu(Azu), azv(Azv),
      cyy(Cyy), cyz(Cyz), czy(Czy), czz(Czz)
  {}
private:
  double ky, kz;
  double ayu, ayv, azu, azv;
  double cyy, cyz, czy, czz;

public:
  void Print( std::ostream &ost ) const;

  friend RKBeamTrajectoryPoint RKBeamtraceOneStep( double, const RKBeamTrajectoryPoint & );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &, 
				const RKBeamFieldIntegral &,
				const RKBeamDeltaFieldIntegral &, 
				const RKBeamDeltaFieldIntegral &, double );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &,
				const RKBeamFieldIntegral & );
};

class RKBeamDeltaFieldIntegral
{
public:
  RKBeamDeltaFieldIntegral( double dKyy, double dKyz, double dKyu,
			    double dKyv, double dKyq,
			    double dKzy, double dKzz, double dKzu,
			    double dKzv, double dKzq )
    : dkyy(dKyy), dkyz(dKyz), dkyu(dKyu), dkyv(dKyv), dkyq(dKyq),
      dkzy(dKzy), dkzz(dKzz), dkzu(dKzu), dkzv(dKzv), dkzq(dKzq)
  {}
private:
  double dkyy, dkyz, dkyu, dkyv, dkyq;
  double dkzy, dkzz, dkzu, dkzv, dkzq;
public:
  void Print( std::ostream &ost ) const;
  
  friend RKBeamTrajectoryPoint RKBeamtraceOneStep( double, const RKBeamTrajectoryPoint & );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &, 
				const RKBeamFieldIntegral &,
				const RKBeamDeltaFieldIntegral &, 
				const RKBeamDeltaFieldIntegral &, double );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &,
				const RKBeamFieldIntegral & );
};

class RKBeamCordParameter
{
public:
  RKBeamCordParameter( ) {}
  RKBeamCordParameter( double X, double Y, double Z,
		       double U, double V, double Q )
    : x(X), y(Y), z(Z), u(U), v(V), q(Q)
  {}
  
  RKBeamCordParameter( const ThreeVector &pos,
		       double U, double V, double Q )
    : x(pos.x()), y(pos.y()), z(pos.z()),
      u(U), v(V), q(Q)
  {}
  
  RKBeamCordParameter( const ThreeVector &pos,
		       const ThreeVector &mom );
private:
  double x, y, z, u, v, q;
public:
  ThreeVector PositionInGlobal( void ) const
  { return ThreeVector( x, y, z ); }
  ThreeVector MomentumInGlobal( void ) const;
  void Print( std::ostream &ost ) const; 
  
  double X( void ) const { return x; }
  double Y( void ) const { return y; }
  double Z( void ) const { return z; }
  double U( void ) const { return u; }
  double V( void ) const { return v; }
  double Q( void ) const { return q; }
  
  friend class RKBeamTrajectoryPoint;
  friend RKBeamTrajectoryPoint 
  RKBeamtraceOneStep( double, const RKBeamTrajectoryPoint & );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &, 
				const RKBeamFieldIntegral &,
				const RKBeamDeltaFieldIntegral &, 
				const RKBeamDeltaFieldIntegral &, double );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &,
				const RKBeamFieldIntegral & );
  friend bool 
  RKBeamcheckCrossing( int, const RKBeamTrajectoryPoint &, 
		       const RKBeamTrajectoryPoint &, RKBeamcalcHitPoint & );
};


class RKBeamTrajectoryPoint
{
public:
  RKBeamTrajectoryPoint( double X, double Y, double Z,
			 double U, double V, double Q,
			 double Dydy, double Dydz, double Dydu, 
			 double Dydv, double Dydq,
			 double Dzdy, double Dzdz, double Dzdu, 
			 double Dzdv, double Dzdq,
			 double Dudy, double Dudz, double Dudu, 
			 double Dudv, double Dudq,
			 double Dvdy, double Dvdz, double Dvdu, 
			 double Dvdv, double Dvdq,
			 double L )
    : r(X,Y,Z,U,V,Q),
      dydy(Dydy), dydz(Dydz), dydu(Dydu), dydv(Dydv), dydq(Dydq),
      dzdy(Dzdy), dzdz(Dzdz), dzdu(Dzdu), dzdv(Dzdv), dzdq(Dzdq),
      dudy(Dudy), dudz(Dudz), dudu(Dudu), dudv(Dudv), dudq(Dudq),
      dvdy(Dvdy), dvdz(Dvdz), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
      l(L)
  {}

  RKBeamTrajectoryPoint( const RKBeamCordParameter &R,
			 double Dydy, double Dydz, double Dydu, 
			 double Dydv, double Dydq,
			 double Dzdy, double Dzdz, double Dzdu, 
			 double Dzdv, double Dzdq,
			 double Dudy, double Dudz, double Dudu, 
			 double Dudv, double Dudq,
			 double Dvdy, double Dvdz, double Dvdu, 
			 double Dvdv, double Dvdq,
			 double L )
    : r(R),
      dydy(Dydy), dydz(Dydz), dydu(Dydu), dydv(Dydv), dydq(Dydq),
      dzdy(Dzdy), dzdz(Dzdz), dzdu(Dzdu), dzdv(Dzdv), dzdq(Dzdq),
      dudy(Dudy), dudz(Dudz), dudu(Dudu), dudv(Dudv), dudq(Dudq),
      dvdy(Dvdy), dvdz(Dvdz), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
      l(L)
  {}
  
  RKBeamTrajectoryPoint( const ThreeVector &pos,
			 double U, double V, double Q,
			 double Dydy, double Dydz, double Dydu, 
			 double Dydv, double Dydq,
			 double Dzdy, double Dzdz, double Dzdu, 
			 double Dzdv, double Dzdq,
			 double Dudy, double Dudz, double Dudu, 
			 double Dudv, double Dudq,
			 double Dvdy, double Dvdz, double Dvdu, 
			 double Dvdv, double Dvdq,
			 double L )
    : r(pos,U,V,Q),
      dydy(Dydy), dydz(Dydz), dydu(Dydu), dydv(Dydv), dydq(Dydq),
      dzdy(Dzdy), dzdz(Dzdz), dzdu(Dzdu), dzdv(Dzdv), dzdq(Dzdq),
      dudy(Dudy), dudz(Dudz), dudu(Dudu), dudv(Dudv), dudq(Dudq),
      dvdy(Dvdy), dvdz(Dvdz), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
      l(L)
  {}
  
  RKBeamTrajectoryPoint( const ThreeVector &pos,
			 const ThreeVector &mom,
			 double Dydy, double Dydz, double Dydu, 
			 double Dydv, double Dydq,
			 double Dzdy, double Dzdz, double Dzdu, 
			 double Dzdv, double Dzdq,
			 double Dudy, double Dudz, double Dudu, 
			 double Dudv, double Dudq,
			 double Dvdy, double Dvdz, double Dvdu, 
			 double Dvdv, double Dvdq,
			 double L )
    : r(pos,mom),
      dydy(Dydy), dydz(Dydz), dydu(Dydu), dydv(Dydv), dydq(Dydq),
      dzdy(Dzdy), dzdz(Dzdz), dzdu(Dzdu), dzdv(Dzdv), dzdq(Dzdq),
      dudy(Dudy), dudz(Dudz), dudu(Dudu), dudv(Dudv), dudq(Dudq),
      dvdy(Dvdy), dvdz(Dvdz), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq),
      l(L)
  {}
  
private:
  RKBeamCordParameter r;
  double dydy, dydz, dydu, dydv, dydq;
  double dzdy, dzdz, dzdu, dzdv, dzdq;
  double dudy, dudz, dudu, dudv, dudq;
  double dvdy, dvdz, dvdu, dvdv, dvdq;
  double l;
public:
  ThreeVector PositionInGlobal( void ) const 
  { return r.PositionInGlobal(); }
  ThreeVector MomentumInGlobal( void ) const
  { return r.MomentumInGlobal(); }
  double PathLength( void ) const { return l; }
  void Print( std::ostream &ost ) const; 
  
  friend RKBeamTrajectoryPoint 
  RKBeamtraceOneStep( double, const RKBeamTrajectoryPoint & );
  friend bool 
  RKBeamcheckCrossing( int, const RKBeamTrajectoryPoint &, 
		       const RKBeamTrajectoryPoint &, RKBeamcalcHitPoint & );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &, 
				const RKBeamFieldIntegral &,
				const RKBeamDeltaFieldIntegral &, 
				const RKBeamDeltaFieldIntegral &, double );
  friend RKBeamDeltaFieldIntegral
  RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &,
				const RKBeamFieldIntegral & );
};

class RKBeamcalcHitPoint {
public:
  RKBeamcalcHitPoint() {}
  RKBeamcalcHitPoint( const ThreeVector &pos, const ThreeVector &mom,
		      double S, double L, 
		      double Dsdy,  double Dsdz,  double Dsdu, 
		      double Dsdv,  double Dsdq, 
		      double Dsdyy, double Dsdyz, double Dsdyu,
		      double Dsdyv, double Dsdyq,
		      double Dsdzy, double Dsdzz, double Dsdzu,
		      double Dsdzv, double Dsdzq,
		      double Dsduy, double Dsduz, double Dsduu,
		      double Dsduv, double Dsduq,
		      double Dsdvy, double Dsdvz, double Dsdvu,
		      double Dsdvv, double Dsdvq,
		      double Dsdqy, double Dsdqz, double Dsdqu,
		      double Dsdqv, double Dsdqq,
		      double Dydy,  double Dydz,  double Dydu, 
		      double Dydv,  double Dydq,
		      double Dzdy,  double Dzdz,  double Dzdu, 
		      double Dzdv,  double Dzdq,
		      double Dudy,  double Dudz,  double Dudu, 
		      double Dudv,  double Dudq,
		      double Dvdy,  double Dvdz,  double Dvdu, 
		      double Dvdv,  double Dvdq )
    : posG(pos), momG(mom), s(S), l(L),
      dsdy(Dsdy), dsdz(Dsdz), dsdu(Dsdu), dsdv(Dsdv), dsdq(Dsdq), 
      dsdyy(Dsdyy), dsdyz(Dsdyz), dsdyu(Dsdyu), dsdyv(Dsdyv), dsdyq(Dsdyq),
      dsdzy(Dsdzy), dsdzz(Dsdzz), dsdzu(Dsdzu), dsdzv(Dsdzv), dsdzq(Dsdzq),
      dsduy(Dsduy), dsduz(Dsduz), dsduu(Dsduu), dsduv(Dsduv), dsduq(Dsduq),
      dsdvy(Dsdvy), dsdvz(Dsdvz), dsdvu(Dsdvu), dsdvv(Dsdvv), dsdvq(Dsdvq),
      dsdqy(Dsdqy), dsdqz(Dsdqz), dsdqu(Dsdqu), dsdqv(Dsdqv), dsdqq(Dsdqq),
      dydy(Dydy), dydz(Dydz), dydu(Dydu), dydv(Dydv), dydq(Dydq), 
      dzdy(Dzdy), dzdz(Dzdz), dzdu(Dzdu), dzdv(Dzdv), dzdq(Dzdq), 
      dudy(Dudy), dudz(Dudz), dudu(Dudu), dudv(Dudv), dudq(Dydq), 
      dvdy(Dvdy), dvdz(Dvdz), dvdu(Dvdu), dvdv(Dvdv), dvdq(Dvdq) 
  {}
  
private:
  ThreeVector posG, momG;
  double s;
  double l;
  double dsdy, dsdz, dsdu, dsdv, dsdq;
  double dsdyy, dsdyz, dsdyu, dsdyv, dsdyq;
  double dsdzy, dsdzz, dsdzu, dsdzv, dsdzq;
  double dsduy, dsduz, dsduu, dsduv, dsduq;
  double dsdvy, dsdvz, dsdvu, dsdvv, dsdvq;
  double dsdqy, dsdqz, dsdqu, dsdqv, dsdqq;
  double dydy, dydz, dydu, dydv, dydq;
  double dzdy, dzdz, dzdu, dzdv, dzdq;
  double dudy, dudz, dudu, dudv, dudq;
  double dvdy, dvdz, dvdu, dvdv, dvdq;
public:
  const ThreeVector &PositionInGlobal( void ) const { return posG; }
  const ThreeVector &MomentumInGlobal( void ) const { return momG; }
  double PositionInLocal( void ) const { return s; } 
  double PathLength( void ) const { return l; }
  double coefY( void ) const { return dsdy; }
  double coefZ( void ) const { return dsdz; }
  double coefU( void ) const { return dsdu; }
  double coefV( void ) const { return dsdv; }
  double coefQ( void ) const { return dsdq; }
  double coefYY( void ) const { return dsdyy; }
  double coefYZ( void ) const { return dsdyz; }
  double coefYU( void ) const { return dsdyu; }
  double coefYV( void ) const { return dsdyv; }
  double coefYQ( void ) const { return dsdyq; }
  double coefZY( void ) const { return dsdzy; }
  double coefZZ( void ) const { return dsdzz; }
  double coefZU( void ) const { return dsdzu; }
  double coefZV( void ) const { return dsdzv; }
  double coefZQ( void ) const { return dsdzq; }
  double coefUY( void ) const { return dsduy; }
  double coefUZ( void ) const { return dsduz; }
  double coefUU( void ) const { return dsduu; }
  double coefUV( void ) const { return dsduv; }
  double coefUQ( void ) const { return dsduq; }
  double coefVY( void ) const { return dsdvy; }
  double coefVZ( void ) const { return dsdvz; }
  double coefVU( void ) const { return dsdvu; }
  double coefVV( void ) const { return dsdvv; }
  double coefVQ( void ) const { return dsdvq; }
  double coefQY( void ) const { return dsdqy; }
  double coefQZ( void ) const { return dsdqz; }
  double coefQU( void ) const { return dsdqu; }
  double coefQV( void ) const { return dsdqv; }
  double coefQQ( void ) const { return dsdqq; }

  double dYdY( void ) const { return dydy; }
  double dYdZ( void ) const { return dydz; }
  double dYdU( void ) const { return dydu; }
  double dYdV( void ) const { return dydv; }
  double dYdQ( void ) const { return dydq; }
  double dZdY( void ) const { return dzdy; }
  double dZdZ( void ) const { return dzdz; }
  double dZdU( void ) const { return dzdu; }
  double dZdV( void ) const { return dzdv; }
  double dZdQ( void ) const { return dzdq; }
  double dUdY( void ) const { return dudy; }
  double dUdZ( void ) const { return dudz; }
  double dUdU( void ) const { return dudu; }
  double dUdV( void ) const { return dudv; }
  double dUdQ( void ) const { return dudq; }
  double dVdY( void ) const { return dvdy; }
  double dVdZ( void ) const { return dvdz; }
  double dVdU( void ) const { return dvdu; }
  double dVdV( void ) const { return dvdv; }
  double dVdQ( void ) const { return dvdq; }

  friend bool 
  RKBeamcheckCrossing( int, const RKBeamTrajectoryPoint &, 
		       const RKBeamTrajectoryPoint &, RKBeamcalcHitPoint & );
  
};

class RKBeamHitPointContainer  
  : public std::vector<std::pair<int,RKBeamcalcHitPoint> >
{
public:
  const RKBeamcalcHitPoint & HitPointOfLayer( int lnum ) const;
  RKBeamcalcHitPoint & HitPointOfLayer( int lnum );
  
  typedef std::vector<std::pair<int,RKBeamcalcHitPoint> >
  ::const_iterator RKBeamHpCIterator;
  
  typedef std::vector<std::pair<int,RKBeamcalcHitPoint> >
  ::iterator RKBeamHpIterator;
  
};

/////////////////////////////////////////////////////////////////////////
inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKBeamFieldIntegral &obj )
{ obj.Print( ost ); return ost; }

inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKBeamDeltaFieldIntegral &obj )
{ obj.Print( ost ); return ost; }

inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKBeamCordParameter &obj ) 
{ obj.Print( ost ); return ost; }
inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKBeamTrajectoryPoint &obj ) 
{ obj.Print( ost ); return ost; }
//////////////////////////////////////////////////////////////////////////
RKBeamFieldIntegral 
RKBeamcalcFieldIntegral( double U, double V, double Q, const ThreeVector &B );

RKBeamFieldIntegral 
RKBeamcalcFieldIntegral( double U, double V, double Q, const ThreeVector &B,
			 const ThreeVector &dBdY,  const ThreeVector &dBdZ );

RKBeamDeltaFieldIntegral
RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &prevPoint,
			      const RKBeamFieldIntegral &intg );

RKBeamDeltaFieldIntegral
RKBeamcalcDeltaFieldIntegral( const RKBeamTrajectoryPoint &prevPoint,
			      const RKBeamFieldIntegral &intg,
			      const RKBeamDeltaFieldIntegral &dIntg1,
			      const RKBeamDeltaFieldIntegral &dIntg2,
			      double StepSize );

RKBeamTrajectoryPoint 
RKBeamtraceOneStep( double StepSize, const RKBeamTrajectoryPoint &prevPoint );


bool RKBeamcheckCrossing( int lnum, const RKBeamTrajectoryPoint &startPoint,
			  const RKBeamTrajectoryPoint &endPoint,
			  RKBeamcalcHitPoint &crossPoint ); 

bool RKBeamtrace( const RKBeamCordParameter &initial,
		  RKBeamHitPointContainer &hitContainer );

bool RKBeamtraceToLast( RKBeamHitPointContainer &hitContainer );

RKBeamHitPointContainer RKBeammakeHPContainer( void );


#endif
