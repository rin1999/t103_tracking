/*
  BH2Hit.hh
*/

#ifndef BH2HIT_H
#define BH2HIT_H

#include "Hodo2Hit.hh"

class BH2Hit : public Hodo2Hit
{
public:
  BH2Hit( HodoRawHit *rhit )
    : Hodo2Hit(rhit), Tofs_(0.)
  {};
  ~BH2Hit(){};
  
public:
  bool calculate( void );
  double Time0( void )  const { return 0.5*(t1_+t2_)+Tofs_; };
  double CTime0( void ) const { return 0.5*(ct1_+ct2_)+Tofs_; };
  
  virtual bool ReCalc( bool applyRecursively=false )
  { return calculate(); };
  
private:
  double Tofs_;
};
#endif
