/*
  TrLTrackHit.cc

  2019/2  K.Shirotori
*/

#include "TrLTrackHit.hh"
#include "TrAnalyzer.hh"

#include <cmath>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

double TrLTrackHit::GetLocalCalPos( void ) const
{
  double angle=Hit_->GetTiltAngle();
  return xcal_*cos(angle*Deg2Rad)+ycal_*sin(angle*Deg2Rad);
} 

bool TrLTrackHit::ReCalc( bool applyRecursively )
{
  // if(applyRecursively)
  //   if(!Hit_->ReCalcTr(applyRecursively))
  //     return false;

  // double pos=GetPosition();
  // double dl=GetDriftLength();

  // if( xl_>pos ) xl_=pos+dl;
  // else         xl_=pos-dl;

  return true;
}
