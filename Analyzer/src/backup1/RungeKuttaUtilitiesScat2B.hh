/*
  RungeKuttaUtilitiesScat2B.hh

  2004/6/28    T.Takahashi

*/

#ifndef RungeKuttaUtilitiesScat2B_h

#define RungeKuttaUtilitiesScat2B_h 1

#include "ThreeVector.hh"

#include <vector>
#include <utility>
#include <iosfwd>

class RKScat2BFieldIntegral;
class RKScat2BDeltaFieldIntegral;
class RKScat2BTrajectoryPoint;
class RKScat2BcalcHitPoint;
class RKScat2BCordParameter;
class RKScat2BHitPointContainer;

class RKScat2BFieldIntegral
{
public:
  RKScat2BFieldIntegral( double Ky, double Kz, 
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

  friend RKScat2BTrajectoryPoint RKScat2BtraceOneStep( double, const RKScat2BTrajectoryPoint & );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &, 
				  const RKScat2BFieldIntegral &,
				  const RKScat2BDeltaFieldIntegral &, 
				  const RKScat2BDeltaFieldIntegral &, double );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &,
				  const RKScat2BFieldIntegral & );
};

class RKScat2BDeltaFieldIntegral
{
public:
  RKScat2BDeltaFieldIntegral( double dKyy, double dKyz, double dKyu,
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
  
  friend RKScat2BTrajectoryPoint RKScat2BtraceOneStep( double, const RKScat2BTrajectoryPoint & );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &, 
				  const RKScat2BFieldIntegral &,
				  const RKScat2BDeltaFieldIntegral &, 
				  const RKScat2BDeltaFieldIntegral &, double );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &,
				  const RKScat2BFieldIntegral & );
};

class RKScat2BCordParameter
{
public:
  RKScat2BCordParameter( ) {}
  RKScat2BCordParameter( double X, double Y, double Z,
			 double U, double V, double Q )
    : x(X), y(Y), z(Z), u(U), v(V), q(Q)
  {}
  
  RKScat2BCordParameter( const ThreeVector &pos,
			 double U, double V, double Q )
    : x(pos.x()), y(pos.y()), z(pos.z()),
      u(U), v(V), q(Q)
  {}
  
  RKScat2BCordParameter( const ThreeVector &pos,
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
  
  friend class RKScat2BTrajectoryPoint;
  friend RKScat2BTrajectoryPoint 
  RKScat2BtraceOneStep( double, const RKScat2BTrajectoryPoint & );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &, 
				  const RKScat2BFieldIntegral &,
				  const RKScat2BDeltaFieldIntegral &, 
				  const RKScat2BDeltaFieldIntegral &, double );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &,
				  const RKScat2BFieldIntegral & );
  friend bool 
  RKScat2BcheckCrossing( int, const RKScat2BTrajectoryPoint &, 
			 const RKScat2BTrajectoryPoint &, RKScat2BcalcHitPoint & );
};


class RKScat2BTrajectoryPoint
{
public:
  RKScat2BTrajectoryPoint( double X, double Y, double Z,
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

  RKScat2BTrajectoryPoint( const RKScat2BCordParameter &R,
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
  
  RKScat2BTrajectoryPoint( const ThreeVector &pos,
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
  
  RKScat2BTrajectoryPoint( const ThreeVector &pos,
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
  RKScat2BCordParameter r;
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
  
  friend RKScat2BTrajectoryPoint 
  RKScat2BtraceOneStep( double, const RKScat2BTrajectoryPoint & );
  friend bool 
  RKScat2BcheckCrossing( int, const RKScat2BTrajectoryPoint &, 
			 const RKScat2BTrajectoryPoint &, RKScat2BcalcHitPoint & );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &, 
				  const RKScat2BFieldIntegral &,
				  const RKScat2BDeltaFieldIntegral &, 
				  const RKScat2BDeltaFieldIntegral &, double );
  friend RKScat2BDeltaFieldIntegral
  RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &,
				  const RKScat2BFieldIntegral & );
};

class RKScat2BcalcHitPoint {
public:
  RKScat2BcalcHitPoint() {}
  RKScat2BcalcHitPoint( const ThreeVector &pos, const ThreeVector &mom,
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
  RKScat2BcheckCrossing( int, const RKScat2BTrajectoryPoint &, 
			 const RKScat2BTrajectoryPoint &, RKScat2BcalcHitPoint & );
  
};

class RKScat2BHitPointContainer  
  : public std::vector<std::pair<int,RKScat2BcalcHitPoint> >
{
public:
  const RKScat2BcalcHitPoint & HitPointOfLayer( int lnum ) const;
  RKScat2BcalcHitPoint & HitPointOfLayer( int lnum );
  
  typedef std::vector<std::pair<int,RKScat2BcalcHitPoint> >
  ::const_iterator RKScat2BHpCIterator;
  
  typedef std::vector<std::pair<int,RKScat2BcalcHitPoint> >
  ::iterator RKScat2BHpIterator;
  
};

/////////////////////////////////////////////////////////////////////////
inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKScat2BFieldIntegral &obj )
{ obj.Print( ost ); return ost; }

inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKScat2BDeltaFieldIntegral &obj )
{ obj.Print( ost ); return ost; }

inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKScat2BCordParameter &obj ) 
{ obj.Print( ost ); return ost; }
inline std::ostream & operator << ( std::ostream &ost, 
                                    const RKScat2BTrajectoryPoint &obj ) 
{ obj.Print( ost ); return ost; }
//////////////////////////////////////////////////////////////////////////
RKScat2BFieldIntegral 
RKScat2BcalcFieldIntegral( double U, double V, double Q, const ThreeVector &B );

RKScat2BFieldIntegral 
RKScat2BcalcFieldIntegral( double U, double V, double Q, const ThreeVector &B,
			  const ThreeVector &dBdY,  const ThreeVector &dBdZ );

RKScat2BDeltaFieldIntegral
RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &prevPoint,
				const RKScat2BFieldIntegral &intg );

RKScat2BDeltaFieldIntegral
RKScat2BcalcDeltaFieldIntegral( const RKScat2BTrajectoryPoint &prevPoint,
				const RKScat2BFieldIntegral &intg,
				const RKScat2BDeltaFieldIntegral &dIntg1,
				const RKScat2BDeltaFieldIntegral &dIntg2,
				double StepSize );

RKScat2BTrajectoryPoint 
RKScat2BtraceOneStep( double StepSize, const RKScat2BTrajectoryPoint &prevPoint );


bool RKScat2BcheckCrossing( int lnum, const RKScat2BTrajectoryPoint &startPoint,
			    const RKScat2BTrajectoryPoint &endPoint,
			    RKScat2BcalcHitPoint &crossPoint ); 

bool RKScat2Btrace( const RKScat2BCordParameter &initial,
		    RKScat2BHitPointContainer &hitContainer );

bool RKScat2BtraceToLast( RKScat2BHitPointContainer &hitContainer );

RKScat2BHitPointContainer RKScat2BmakeHPContainer( void );

#endif
