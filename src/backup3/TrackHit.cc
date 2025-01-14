/*
  TrackHit.cc
*/

#include "TrackHit.hh"
#include "DCGeomMan.hh"
//#include "SksObjectId.hh"
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


////////////
// for DS //
////////////

// std::size_t TrackHit::DSSave( unsigned int *bufp ) const
// {
//   DSFormat data;
//   std::size_t size = sizeof(data)/sizeof(int);
//   if( sizeof(data)%sizeof(int) ) size+=1;

//   data.header_=DSMakeHeader( DSMakeHeaderID( SksObjTrackHit ), size );
//   data.pkey_=dchitp_->DSGetSaveKey();
//   data.gxcal_=calGPos_.x(); data.gycal_=calGPos_.y(); data.gzcal_=calGPos_.z();
//   data.callpos_=calLPos_;

//   std::memcpy( bufp, &data, sizeof(data) );
//   return size;
// }

// bool TrackHit::DSRestore( unsigned int *bufp,
// 			  DCLocalTrack *trIn, DCLocalTrack *trOut )
// {
//   static const std::string funcname = "[TrackHit::DSRestore]";
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
  
//   calGPos_.setX( data.gxcal_ );  calGPos_.setY( data.gycal_ );
//   calGPos_.setZ( data.gzcal_ );  calLPos_=data.callpos_;

//   DCLTrackHit *hit=0;
//   switch( HCId ){
//   case SksObjLocalTrackBIC:
//   case SksObjLocalTrackSIC:
//     hit=trIn->GetHit(seqNum); break;
//   case SksObjLocalTrackBOC:
//   case SksObjLocalTrackSOC:
//     hit=trOut->GetHit(seqNum); break;
//   default:
//     {
//       std::ostringstream mess;
//       mess << funcname << ": invalid HCId" << HCId;
//       std::cerr << mess.str() << std::endl;
//       throw std::invalid_argument(mess.str());
//     }
//   }
//   if(hit){
//     dchitp_=hit;
//     return true;
//   }
//   else
//     return false;

// }
