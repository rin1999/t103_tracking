/*
 RawData.cc
*/

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include "DetectorID.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "Delete.hh"

#ifdef MemoryLeak
debug::Counter RawData::sm_counter("RawData");
#endif

enum EDCDataType
  {
    kDCDataLeading,
    kDCDataTrailing,
    kNDCDataType
  };

RawData::RawData():     
  GCRawHC(0),
  BH1RawHC(0),
  BH2RawHC(0),
  BACRawHC(0),
  TGTRawHC(0),
  TOFRawHC(0),
  LCRawHC(0),
  ACRawHC(),
  SP0RawHC(),

  BcInRawHC(),
  BcOutRawHC(),
  SdcInRawHC(),
  SdcOutRawHC(),

  MiscRawHC(0),
  MatrixRawHC(0)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

RawData::~RawData()
{
  clearAll();
#ifdef MemoryLeak
  --sm_counter;
#endif
}

bool RawData::AddHodoRawHit(HodoRHitContainer& cont,
			    int DetId,
			    int Plane,
			    int Seg,
			    int AorT,
			    int UorD,
			    int Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp(Data);
      else          p->SetAdcDown(Data);
    }else{
      if( UorD==0 ) p->SetTdcUp(Data);
      else          p->SetTdcDown(Data);
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddDCRawHit(DCRHitContainer& cont,
			  int Plane,
			  int Wire,
			  int Tdc,
			  int type)
{
  static const std::string funcname = "[RawData::AddDCRawHit]";

  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->PlaneId()==Plane &&
	q->WireId()==Wire ){
       p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Plane, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
   switch(type)
      {
      case kDCDataLeading:
	p->SetTdc( Tdc );
	break;
      case kDCDataTrailing:
	p->SetTrailing(Tdc);
	break;
      default:
	break;
      }
   return true;
  }else{
    std::cerr << funcname << ": new fail. PlaneId="
              << Plane << " WireId="
              << Wire << std::endl;
    return false;
  }
}

void
RawData::clearAll()
{
  std::for_each(GCRawHC.begin(), GCRawHC.end(), hddaq::Delete());
  GCRawHC.clear();

  std::for_each(BH1RawHC.begin(), BH1RawHC.end(), hddaq::Delete());
  BH1RawHC.clear();

  std::for_each(BH2RawHC.begin(), BH2RawHC.end(), hddaq::Delete());
  BH2RawHC.clear();

  std::for_each(BACRawHC.begin(), BACRawHC.end(), hddaq::Delete());
  BACRawHC.clear();

  std::for_each(TGTRawHC.begin(), TGTRawHC.end(), hddaq::Delete());
  TGTRawHC.clear();

  std::for_each(TOFRawHC.begin(), TOFRawHC.end(), hddaq::Delete());
  TOFRawHC.clear();

  std::for_each(LCRawHC.begin(), LCRawHC.end(), hddaq::Delete());
  LCRawHC.clear();

  for( int l=0; l<=NumOfLayersAc; ++l ){
    for_each( ACRawHC[l].begin(),  ACRawHC[l].end(),  hddaq::Delete());
    ACRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSP0; ++l ){
    for_each( SP0RawHC[l].begin(),  SP0RawHC[l].end(),  hddaq::Delete());
    SP0RawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcIn; ++l){
    for_each( BcInRawHC[l].begin(),  BcInRawHC[l].end(),   hddaq::Delete());
    BcInRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcOut; ++l){
    for_each( BcOutRawHC[l].begin(), BcOutRawHC[l].end(),  hddaq::Delete());
    BcOutRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcIn; ++l){
    for_each( SdcInRawHC[l].begin(),  SdcInRawHC[l].end(),   hddaq::Delete());
    SdcInRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcOut; ++l){
    for_each( SdcOutRawHC[l].begin(), SdcOutRawHC[l].end(),  hddaq::Delete());
    SdcOutRawHC[l].clear();
  }

  std::for_each(MiscRawHC.begin(), MiscRawHC.end(), hddaq::Delete());
  MiscRawHC.clear();

  std::for_each(MatrixRawHC.begin(), MatrixRawHC.end(), hddaq::Delete());
  MatrixRawHC.clear();
    
  return;
}

bool
RawData::DecodeHits()
{
  hddaq::unpacker::UnpackerManager& gUnpacker
    = hddaq::unpacker::GUnpacker::get_instance();
  clearAll();

  bool MHTDCflag=false;
  MHTDCflag=ConfMan::GetConfManager()->GetMHTDCFlag();

  int offset=0;
  if( MHTDCflag ) offset = 64000;//(0x10000 - 0x600)
  else offset = 0;

  //std::cout << "Offset=" << offset << std::endl; 

  //GC
  for( int AorT=0; AorT<2;  ++AorT ){
    int nhit = gUnpacker.get_entries( DetIdGC, 0, 0, 0, AorT );
    if( nhit>0 ){
      int data = gUnpacker.get( DetIdGC, 0, 0, 0, AorT );
      AddHodoRawHit( GCRawHC, DetIdGC, 0, 0, AorT, 0, data );
    }
    else continue;
  }  
  
  //BH1
  for( int seg=0; seg<NumOfSegBH1; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdBH1, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdBH1, 0, seg, UorD, AorT );
	  AddHodoRawHit( BH1RawHC, DetIdBH1, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }

  //BH2 
  for( int seg=0; seg<NumOfSegBH2; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdBH2, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdBH2, 0, seg, UorD, AorT );
	  AddHodoRawHit( BH2RawHC, DetIdBH2, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }

  //BAC
  for( int seg=0; seg<NumOfSegBAC; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdBAC, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdBAC, 0, seg, 0, AorT );
	AddHodoRawHit( BACRawHC, DetIdBAC, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }

  //TGT
  for( int seg=0; seg<NumOfSegTGT; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdTGT, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdTGT, 0, seg, 0, AorT );
	AddHodoRawHit( TGTRawHC, DetIdTGT, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }

  //TOF
  for( int seg=0; seg<NumOfSegTOF; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdTOF, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdTOF, 0, seg, UorD, AorT );
	  AddHodoRawHit( TOFRawHC, DetIdTOF, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }

  //LC
  for( int seg=0; seg<NumOfSegLC; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdLC, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdLC, 0, seg, UorD, AorT );
	  AddHodoRawHit( LCRawHC, DetIdLC, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }
  
  //AC
  for( int seg=0; seg<NumOfSegAC; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdAC, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdAC, 0, seg, UorD, AorT );
	  AddHodoRawHit( ACRawHC[UorD+1], DetIdAC, UorD+1, seg, AorT, 0, data );
	}
	else continue;
      }
    }
  }
  
  //SP0
//   for(int plane=0; plane<NumOfLayersSP0; ++plane ){
//     for( int seg=0; seg<NumOfSegSP0; ++seg ){
//       for( int AorT=0; AorT<2;  ++AorT ){
// 	for( int UorD=0; UorD<2; ++UorD ){
// 	  int nhit = gUnpacker.get_entries( DetIdSP0, plane, seg, UorD, AorT );
// 	  if( nhit>0 ){
// 	    int data = gUnpacker.get( DetIdAC, 0, seg, UorD, AorT );
// 	    AddHodoRawHit( SP0RawHC[plane+1], DetIdSP0, plane+1, seg, AorT, 0, data );
// 	  }
// 	  else continue;
// 	}
//       }
//     }
//   }

  // BC1&BC2 MWPC
  for(int plane=0; plane<NumOfLayersBcIn; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC1; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC1, plane, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int leading = gUnpacker.get( DetIdBC1, plane, 0, wire, 0, i );
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 leading, kDCDataLeading );
	    int trailing = gUnpacker.get( DetIdBC1, plane, 0, wire, 1, i );
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 trailing, kDCDataTrailing );
	  }
	}
	else continue; 
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC2; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int leading = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0, i );
	    int trailing = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 1, i );
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 leading, kDCDataLeading); 
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 trailing, kDCDataTrailing); 
	  }
	}	
	else continue;
      }
    }
  }
  
  // BC3&BC4 MWDC
  for(int plane=0; plane<NumOfLayersBcOut; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC3; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC3, plane, 0, wire, 2 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = ( gUnpacker.get( DetIdBC3, plane, 0, wire, 2, i ) & 0xffff);
	    AddDCRawHit( BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data-offset );
	  }
	}
	else continue; 
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC4; ++wire){
	int nhit = gUnpacker.get_entries( DetIdK6BDC, plane-NumOfLayersBc, 0, wire, 2 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = ( gUnpacker.get( DetIdK6BDC, plane-NumOfLayersBc, 0, wire, 2, i ) & 0xffff);
	    AddDCRawHit( BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data-offset );
	  }
	}	
	else continue;
      }
    }
  }

  // SDC1&SDC2 MWDC
  for(int plane=0; plane<NumOfLayersSdcIn; ++plane ){
    if( plane<NumOfLayersSdc-2 ){
      for(int wire=0; wire<MaxWireSDC1; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC1, plane, 0, wire, 2 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = (gUnpacker.get( DetIdSDC1, plane, 0, wire, 2, i ) & 0xffff);
	    AddDCRawHit( SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data-offset );
	  }
	}
	else continue; 
      }
    }
    else{
      for(int wire=0; wire<MaxWireSDC2; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC2, plane-(NumOfLayersSdc-2), 0, wire, 2 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = (gUnpacker.get( DetIdSDC2, plane-(NumOfLayersSdc-2), 0, wire, 2, i ) & 0xffff);
	    AddDCRawHit( SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data-offset );
	  }
	}	
	else continue;
      }
    }
  }

  // SDC3&SDC4
  for(int plane=0; plane<NumOfLayersSdcOut; ++plane ){
    if( plane<NumOfLayersSdc ){
      int NumOfWireSDC3;
      if( plane==1 || plane==4 ) NumOfWireSDC3 =108;
      else NumOfWireSDC3 =120;
      for(int wire=0; wire<NumOfWireSDC3; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC3, plane, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC3, plane, 0, wire, 0, i );
 	    AddDCRawHit( SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1, data );
	  }
	}
	else continue;
      }
    }
    else{
      int NumOfWireSDC4;
      if( plane==(NumOfLayersSdc+1) || plane==(NumOfLayersSdc+4) ) NumOfWireSDC4 =108;
      else NumOfWireSDC4 =120;
      for(int wire=0; wire<NumOfWireSDC4; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC4, plane-NumOfLayersSdc, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC4, plane-NumOfLayersSdc, 0, wire, 0 ,i );
	    AddDCRawHit( SdcOutRawHC[plane+1],  plane+PlMinSdcOut, NumOfWireSDC4-wire, data );
	    //	    AddDCRawHit( SdcOutRawHC[plane+1],  plane+PlMinSdcOut, wire+1, data );
	  }
	}
	else continue; 
      }
    }
  }

  //Misc
  for( int seg=0; seg<NumOfMisc; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdMisc, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdMisc, 0, seg, 0, AorT );
	AddHodoRawHit( MiscRawHC, DetIdMisc, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }

  //Matrix
  for( int seg=0; seg<NumOfMatrix; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdMatrix, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdMatrix, 0, seg, 0, AorT );
	AddHodoRawHit( MatrixRawHC, DetIdMatrix, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }


  return true;
}

const HodoRHitContainer& RawData::GetGCRawHC() const
{
  return GCRawHC;
}

const HodoRHitContainer& RawData::GetBH1RawHC() const
{
  return BH1RawHC;
}

const HodoRHitContainer& RawData::GetBH2RawHC() const
{
  return BH2RawHC;
}

const HodoRHitContainer& RawData::GetBACRawHC() const
{
  return BACRawHC;
}

const HodoRHitContainer& RawData::GetTGTRawHC() const
{
  return TGTRawHC;
}

const HodoRHitContainer& RawData::GetTOFRawHC() const
{
  return TOFRawHC;
}

const HodoRHitContainer& RawData::GetLCRawHC() const
{
  return LCRawHC;
}

const HodoRHitContainer& RawData::GetACRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersAc ) layer=0;
  return ACRawHC[layer];
}

const HodoRHitContainer& RawData::GetSP0RawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSP0 ) layer=0;
  return SP0RawHC[layer];
}

const DCRHitContainer & RawData::GetBcInRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return BcInRawHC[layer];
}

const DCRHitContainer & RawData::GetBcOutRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcOut ) layer=0;
  return BcOutRawHC[layer];
}

const DCRHitContainer & RawData::GetSdcInRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcIn ) layer=0;
  return SdcInRawHC[layer];
}

const DCRHitContainer & RawData::GetSdcOutRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcOut ) layer=0;
  return SdcOutRawHC[layer];
}

const HodoRHitContainer& RawData::GetMiscRawHC() const
{
  return MiscRawHC;
}

const HodoRHitContainer& RawData::GetMatrixRawHC() const
{
  return MatrixRawHC;
}

