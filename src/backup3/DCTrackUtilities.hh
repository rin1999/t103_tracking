/*
  DCTrackUtilities.hh

  2005/6/28
*/

#ifndef DCTrackUtilities_h

#define DCTrackUtilities_h

#include "ThreeVector.hh"

class DCLocalTrack;


bool HasHitsAlongStraightLine( DCLocalTrack * TrSdcIn, 
			       const ThreeVector & bPos,
			       const ThreeVector & bMom );

#endif
