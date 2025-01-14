/*
  PrimInfo.hh
  
  2016/2  K.Shirotori
*/

#ifndef PrimInfo_h 
#define PrimInfo_h

#include <cstddef>
#include <vector>

class PrimInfo
{
private:
  double x_, y_, z_;
  double pb_, ub_, vb_;
  double abmom_;
  double m1_, p1_, theta1_, phi1_, thetaCM1_, phiCM1_;
  double m2_, p2_, theta2_, phi2_, thetaCM2_, phiCM2_;

public:
  PrimInfo()
    : x_(0.0), y_(0.0), z_(0.0),
      pb_(0.0), ub_(0.0), vb_(0.0),
      abmom_(0.0),
      m1_(0.0), p1_(0.0), 
      theta1_(0.0), phi1_(0.0), thetaCM1_(0.0), phiCM1_(0.0),
      m2_(0.0), p2_(0.0), 
      theta2_(0.0), phi2_(0.0), thetaCM2_(0.0), phiCM2_(0.0)
  {};
  ~PrimInfo() {};

public:
  void SetVertX( double x ) { x_=x; }
  void SetVertY( double y ) { y_=y; }
  void SetVertZ( double z ) { z_=z; }
  void SetBeamMom( double pb ) { pb_=pb; }
  void SetBeamU( double ub ) { ub_=ub; }
  void SetBeamV( double vb ) { vb_=vb; }
  void SetAnaBeamMom( double abmom ) { abmom_=abmom; }
  void SetMass1( double m1 ) { m1_=m1; }
  void SetMom1( double p1 ) { p1_=p1; }
  void SetTheta1( double theta1 ) { theta1_=theta1; }
  void SetPhi1( double phi1 ) { phi1_=phi1; }
  void SetThetaCM1( double thetaCM1 ) { thetaCM1_=thetaCM1; }
  void SetPhiCM1( double phiCM1 ) { phiCM1_=phiCM1; }
  void SetMass2( double m2 ) { m2_=m2; }
  void SetMom2( double p2 ) { p2_=p2; }
  void SetTheta2( double theta2 ) { theta2_=theta2; }
  void SetPhi2( double phi2 ) { phi2_=phi2; }
  void SetThetaCM2( double thetaCM2 ) { thetaCM2_=thetaCM2; }
  void SetPhiCM2( double phiCM2 ) { phiCM2_=phiCM2; }

  double GetVertX( void ) const { return x_; };
  double GetVertY( void ) const { return y_; };
  double GetVertZ( void ) const { return z_; };
  double GetBeamMom( void ) const { return pb_; };
  double GetBeamU( void ) const { return ub_; };
  double GetBeamV( void ) const { return vb_; };
  double GetAnaBeamMom( void ) const { return abmom_; };
  double GetMass1( void ) const { return m1_; };
  double GetMom1( void ) const { return p1_; };
  double GetTheta1( void ) const { return theta1_; };
  double GetPhi1( void ) const { return phi1_; };
  double GetThetaCM1( void ) const { return thetaCM1_; };
  double GetPhiCM1( void ) const { return phiCM1_; };
  double GetMass2( void ) const { return m2_; };
  double GetMom2( void ) const { return p2_; };
  double GetTheta2( void ) const { return theta2_; };
  double GetPhi2( void ) const { return phi2_; };
  double GetThetaCM2( void ) const { return thetaCM2_; };
  double GetPhiCM2( void ) const { return phiCM2_; };
};
#endif
