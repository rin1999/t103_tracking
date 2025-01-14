/*
  TrackHit.cc

  2018/12  K.Shirotori
*/

#include "TrackHit.hh"
#include "DCGeomMan.hh"
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
