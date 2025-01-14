/*
  DCLTrackHit.cc

  2018/12  K.Shirotori
*/

#include "DCLTrackHit.hh"
#include "DCAnalyzer.hh"

#include <cmath>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

double DCLTrackHit::GetLocalCalPos( void ) const
{
  double angle=Hit_->GetTiltAngle();
  return xcal_*cos(angle*Deg2Rad)+ycal_*sin(angle*Deg2Rad);
} 

bool DCLTrackHit::ReCalc( bool applyRecursively )
{
  if(applyRecursively)
    if(!Hit_->ReCalcDC(applyRecursively))
      return false;

  double wp=GetWirePosition();
  double dl=GetDriftLength();

  if( xl_>wp ) xl_=wp+dl;
  else         xl_=wp-dl;

  return true;
}
