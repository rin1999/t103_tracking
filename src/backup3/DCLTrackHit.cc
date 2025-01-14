/*
  DCLTrackHit.cc
*/

#include "DCLTrackHit.hh"
//#include "SksObjectId.hh"
#include "DCAnalyzer.hh"

#include <cmath>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

#ifdef MemoryLeak
debug::Counter DCLTrackHit::sm_counter("DCLTrackHit");
#endif

double DCLTrackHit::GetLocalCalPos( void ) const
{
  double angle=Hit_->GetTiltAngle();
  return xcal_*cos(angle*Deg2Rad)+ycal_*sin(angle*Deg2Rad);
} 

bool DCLTrackHit::ReCalc( bool applyRecursively )
{
  if(applyRecursively)
    if(!Hit_->ReCalcDC(applyRecursively)
       || !Hit_->ReCalcMWPC(applyRecursively))
      return false;

  double wp=GetWirePosition();
  double dl=GetDriftLength();

  if( xl_>wp ) xl_=wp+dl;
  else         xl_=wp-dl;

  return true;
}

////////////
// for DS //
////////////

// std::size_t DCLTrackHit::DSSave( unsigned int *bufp ) const
// {
//   DSFormat data;
//   std::size_t size = sizeof(data)/sizeof(int);
//   if( sizeof(data)%sizeof(int) ) size+=1;

//   data.header_=DSMakeHeader( DSMakeHeaderID( SksObjDCLTrackHit ), size );
//   data.pkey_=Hit_->DSGetSaveKey();
//   data.xl_=xl_; data.xcal_=xcal_; data.ycal_=ycal_;

//   std::memcpy( bufp, &data, sizeof(data) );
//   return size;
// }

// bool DCLTrackHit::DSRestore( unsigned int *bufp, DCAnalyzer *DCana )
// {
//   static const std::string funcname = "[DCLTrackHit::DSRestore]";
//   DSFormat data;
//   bool ret=false;

//   std::memcpy( &data, bufp, sizeof(data) );

//   if( DSGetSize(data.header_)==0 ){
//     std::ostringstream mess;
//     mess << funcname << ": Size=0";
//     std::cerr << mess.str() << std::endl;
//     throw std::invalid_argument(mess.str());
//   }

//   int header = DSGetHeaderID( data.pkey_ );
//   int seqNum = DSGetSeqNum( data.pkey_ );
//   int HCId   = DSGetObjectID( header );

//   xl_=data.xl_; xcal_=data.xcal_; ycal_=data.ycal_;

//   DCHit *hit=0;
//   int layer=DSGetObjectSubID( header );

// #if 0
//   std::cout << funcname << ": HCId=" << HCId
// 	    << " Layer=" << layer << " SeqNum="
// 	    << seqNum << std::endl;
// #endif

//   switch( HCId ){
//   case SksObjBdcInHC:
// #if 0
//     std::cout << funcname << ": SksObjBdcInHC " <<
//       DCana->GetBdcInHC(layer).size() << std::endl;
// #endif
//     hit=(DCana->GetBdcInHC(layer))[seqNum]; break;
//   case SksObjBdcOutHC:
//     hit=(DCana->GetBdcOutHC(layer))[seqNum]; break;
//   case SksObjSdcInHC:
//     hit=(DCana->GetSdcInHC(layer))[seqNum]; break;
//   case SksObjSdcOutHC:
//     hit=(DCana->GetSdcOutHC(layer))[seqNum]; break;
//   default:
//     {
//       std::ostringstream mess;
//       mess << funcname << ": invalid HCId " << HCId;
//       std::cerr << mess.str() << std::endl;
//       throw std::invalid_argument(mess.str());
//     }
//   }
// #if 0
//   std::cout << funcname << ":::::" << std::endl;
// #endif

//   if(hit){
//     Hit_=hit; Hit_->RegisterHits(this);
//     return true;
//   }
//   else
//     return false;
// }

