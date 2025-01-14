/*
  TrackHit.cc
*/

#include "TrackHit.hh"
#include "TrGeomMan.hh"
#include "TrLocalTrack.hh"

#include <cstring>
#include <stdexcept>
#include <sstream>

TrackHit::TrackHit( TrLTrackHit *hit )
  : trhitp_(hit)
{

}

TrackHit::~TrackHit()
{
} 
