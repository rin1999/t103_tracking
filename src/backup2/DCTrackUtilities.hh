/*
  DCTrackUtilities.hh

  2018/12  K.Shirotori
*/

#ifndef DCTrackUtilities_h

#define DCTrackUtilities_h

#include "ThreeVector.hh"

class DCLocalTrack;


bool HasHitsAlongStraightLine( DCLocalTrack * TrSdcIn, 
			       const ThreeVector & bPos,
			       const ThreeVector & bMom );

#endif
