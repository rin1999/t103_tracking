/*
  Hodo1Hit.cc
*/

#include "HodoHit.hh"

#include "ConfMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <stdexcept>
#include <sstream>


HodoHit::HodoHit( HodoRawHit *rhit )
  : raw_(rhit), Status_(false)
{}

HodoHit::~HodoHit()
{}

bool HodoHit::calculate( void )
{
  static const std::string funcname = "[HodoHit::calculate]";
  
  Status_=false;

  return Status_=true;
}
