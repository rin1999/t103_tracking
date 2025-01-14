/*
  TrTrackUtilities.hh

  2019/2  K.Shirotori
*/

#ifndef TrTrackUtilities_h

#define TrTrackUtilities_h

#include "ThreeVector.hh"

class DCLocalTrack;


bool HasHitsAlongStraightLine( DCLocalTrack * TrSdcIn, 
			       const ThreeVector & bPos,
			       const ThreeVector & bMom );

#endif
