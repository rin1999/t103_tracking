/*
  TrackHit.cc

  2024/04  K.Shirotori
*/

#include "TrackHit.hh"
#include "GeomMan.hh"
#include "DCLocalTrack.hh"

#include <cstring>
#include <stdexcept>
#include <sstream>

TrackHit::TrackHit( DCLTrackHit *hit )
  : dchitp_(hit)
{

}

TrackHit::~TrackHit()
{
} 
