/*
  RawData.cc

  2018/10  K.Shirotori
*/

#include "RawData.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "CLHEP/Random/Randomize.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "Decoder.hh"
#include "Common_val_drs4.hh"
#include "Datadrs.hh"
#include "DetectorID.hh"
#include "TrRawHit.hh"
#include "HodoRawHit.hh"
#include "HRTdcRawHit.hh"
#include "TemplateLib.hh"

#include "ConfMan.hh"
#include "CMapMan.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
  SFTRHC(), T0RHC(0), HRTRHC(0)
{}

RawData::~RawData()
{
  clearAll();
}

bool RawData::AddTrRHit( TrRHitContainer& cont,
			 int Layer, int Fiber, int UorD, 
			 int AdcH, int AdcL,
			 IntVec lTdc, IntVec tTdc )
{
  static const std::string funcname = "[RawData::AddTrRHit]";
  
  TrRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    TrRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
  	q->FiberId()==Fiber &&
  	q->GetUorD()==UorD ){
      p=q; break;
    }
  }
  if(!p){
    p = new TrRawHit( Layer, Fiber, UorD );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetAdcH( AdcH );
    p->SetAdcL( AdcL );
    p->SetlTdc( lTdc );
    p->SettTdc( tTdc );

    // std::cout << "***********************" << std::endl;
    // std::cout<< Fiber << "::" << AdcH << std::endl;
    // std::cout<< Fiber << "::" << AdcL << std::endl;

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int Layer, int Seg, int UorD, 
			   int EventNum,
			   DoubleVec Waveform, 
			   DoubleVec Adc, DoubleVec Amp, 
			   DoubleVec Bl, DoubleVec PeakX,
			   uIntVec Tdc, uIntVec Dt, 
			   uIntVec Tdc_2nd, uIntVec Dt_2nd,
			   uIntVec Width, 
			   unsigned int L1_Tdc, unsigned int L1_Tdc1 )
{
  static const std::string funcname = "[RawData::AddHodoRHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->LayerId()==Layer &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Layer, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetEvNum( EventNum );
    p->SetUorD( UorD );
    p->SetWaveform( Waveform );
    p->SetAdc( Adc );
    p->SetAmp( Amp );
    p->SetBl( Bl );
    p->SetPeakX( PeakX );
    p->SetTdc( Tdc );
    p->SetDt( Dt );
    p->SetTdc_2nd( Tdc_2nd );
    p->SetDt_2nd( Dt_2nd );
    p->SetWidth( Width );
    p->SetL1_Tdc( L1_Tdc );
    p->SetL1_Tdc1( L1_Tdc1 );

    //std::cout<< "Ch=" << Seg << " :TDC=" << Tdc.size() << std::endl;

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddHRTdcRHit( HRTdcRHitContainer& cont,
			    int DetId, int Layer, int Seg, int UorD, 
			    int EventNum,
			    uIntVec lTdc, uIntVec tTdc )
{
  static const std::string funcname = "[RawData::AddHRTdcRHit]";
  
  HRTdcRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HRTdcRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->LayerId()==Layer &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HRTdcRawHit( DetId, Layer, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetEvNum( EventNum );
    p->SetUorD( UorD );
    p->SetlTdc( lTdc );
    p->SettTdc( tTdc );

    //std::cout<< "Ch=" << Seg << " :TDC=" << lTdc.size() << std::endl;

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

void RawData::clearAll()
{ 

  for( int l=0; l<=PlMaxSFT; ++l){
    for_each( SFTRHC[l].begin(),  SFTRHC[l].end(), DeleteObject());
    SFTRHC[l].clear();
  }
  
  std::for_each(T0RHC.begin(), T0RHC.end(), DeleteObject());
  T0RHC.clear();

  std::for_each(HRTRHC.begin(), HRTRHC.end(), DeleteObject());
  HRTRHC.clear();

  return;
}

bool RawData::DecodeRawHits( Decoder& gDec )
{
  clearAll();

  ConfMan *confMan = ConfMan::GetConfManager();
  CMapMan *mapMan=ConfMan::GetConfManager()->GetCMapManager();

  /////////////////////////////////////////////////////////
  //Fiber decoding
  int FiberAdcHigh[NumOfFiber], FiberAdcLow[NumOfFiber];
  std::vector<std::vector<int> > FiberTdcLeading;
  std::vector<std::vector<int> > FiberTdcTrailing;

  FiberTdcLeading.resize(NumOfFiber);
  FiberTdcTrailing.resize(NumOfFiber);

  gDec.GetFiberData( FiberAdcHigh,
		     FiberAdcLow,
   		     FiberTdcLeading,
   		     FiberTdcTrailing );

  for( int ich=0; ich<NumOfFiber; ich++ ){
    int Caddr = 1, Naddr = 1, Aaddr = -1;
    Aaddr = ich;
    int DetId, PlId, SegId, UorD;
    bool status;

    // std::cout << "***********************" << std::endl;
    // std::cout<< ich << "::" << FiberAdcHigh[ich] << std::endl;
    // std::cout<< ich << "::" << FiberAdcLow[ich] << std::endl;

    status = mapMan->GetLogical(Caddr,Naddr,Aaddr,DetId,PlId,SegId,UorD);
    if(status){
      AddTrRHit( SFTRHC[PlId+1], PlId, SegId, UorD, 
      		 FiberAdcHigh[ich], FiberAdcLow[ich], 
      		 FiberTdcLeading[ich], FiberTdcTrailing[ich] );

#if 0
      std::cout << "***********************" << std::endl;
      std::cout << "Caddr=" << Caddr
		<< " Naddr=" << Naddr
		<< " Aaddr=" << Aaddr
		<< " --> "
		<< "DetId=" << DetId 
		<< " PlId=" << PlId 
		<< " SegId=" << SegId
		<< " UorD=" << UorD 
		<< std::endl;
#endif
    }
  }

  /////////////////////////////////////////////////////////
  //T0 decoding (DRS4)
  DRS::dataDrs sDRSData;
  gDec.decode(sDRSData);
  gDec.decodeADC(sDRSData);
  gDec.decodeTDC(sDRSData);

  //std::cout<< "DRS4 Event Num = " << sDRSData.eventNum << std::endl;

  for(int ich = 0; ich<DRS::NofChModule*DRS::NofBoards; ++ich){
    int Caddr = 1, Naddr = 2, Aaddr = -1;
    Aaddr = ich;
    int DetId, PlId, SegId, UorD;
    bool status;
    
    status = mapMan->GetLogical(Caddr,Naddr,Aaddr,DetId,PlId,SegId,UorD);
    if(status){
      AddHodoRHit( T0RHC, DetId, PlId, SegId, UorD, 
		   sDRSData.eventNum,
		   sDRSData.data_wf[ich],
		   sDRSData.data_adc[ich], sDRSData.data_amp[ich], 
		   sDRSData.data_bl[ich], sDRSData.data_peakx[ich], 
		   sDRSData.data_tdc[ich], sDRSData.data_dt[ich], 
		   sDRSData.data_tdc_2nd[ich], sDRSData.data_dt_2nd[ich], 
		   sDRSData.data_width[ich], 
		   sDRSData.l1_tdc[0], sDRSData.l1_tdc1[0] );
    }
  
#if 0
      std::cout << "***********************" << std::endl;
      std::cout << "Caddr=" << Caddr
		<< " Naddr=" << Naddr
		<< " Aaddr=" << Aaddr
		<< " --> "
		<< "DetId=" << DetId 
		<< " PlId=" << PlId 
		<< " SegId=" << SegId
		<< " UorD=" << UorD 
		<< std::endl;
#endif
  }

  /////////////////////////////////////////////////////////
  //T0 decoding (HR-TDC)
  MZN_HRTDC::Data_t MZNData;
  MZNData = gDec.GetMznHRTdcData();

  std::vector<std::vector<unsigned int> > HRTdcLeading;
  std::vector<std::vector<unsigned int> > HRTdcTrailing;
  HRTdcLeading.resize(NumOfHRTDC);
  HRTdcTrailing.resize(NumOfHRTDC);

  //std::cout<< "MZN Event Num = " << MZNData.EventNum << std::endl;

  if( MZNData.Nhit>0 || MZNData.tNhit>0 ){
    for( int ich=0; ich<NumOfHRTDC; ich++ ){
      int Caddr = 1, Naddr = 3, Aaddr = -1;
      Aaddr = ich;
      int DetId, PlId, SegId, UorD;
      bool status;

      for( int j=0; j<MZNData.Nhit; j++ ){
	if( ich==MZNData.Ch[j] ){
	  HRTdcLeading[ich].push_back(MZNData.tdc_leading[j]);
	}
      }
      for( int j=0; j<MZNData.tNhit; j++ ){
	if( ich==MZNData.tCh[j] ){
	  HRTdcTrailing[ich].push_back(MZNData.tdc_trailing[j]);
	}
      }
      
      status = mapMan->GetLogical(Caddr,Naddr,Aaddr,DetId,PlId,SegId,UorD);
      
      if(status){
	AddHRTdcRHit( HRTRHC, DetId, PlId, SegId, UorD, 
		      MZNData.EventNum,
		      HRTdcLeading[ich], HRTdcTrailing[ich]);
#if 0
      std::cout << "***********************" << std::endl;
      std::cout << "Caddr=" << Caddr
		<< " Naddr=" << Naddr
		<< " Aaddr=" << Aaddr
		<< " --> "
		<< "DetId=" << DetId 
		<< " PlId=" << PlId 
		<< " SegId=" << SegId
		<< " UorD=" << UorD 
		<< std::endl;
#endif
      }
    }
  }

  return true;
}

const TrRHitContainer & RawData::GetSFTRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxSFT ) layer=0;
  return SFTRHC[layer];
}

const HodoRHitContainer& RawData::GetT0RHC() const
{
  return T0RHC;
}

const HRTdcRHitContainer& RawData::GetHRTRHC() const
{
  return HRTRHC;
}

unsigned int RawData::getBigEndian32(const char* b = NULL){
    //std::cout << "size of b " << sizeof(b) << std::endl;
    return ((b[0] << 24) & 0xff000000) |
           ((b[1] << 16) & 0x00ff0000) |
           ((b[2] <<  8) & 0x0000ff00) |
           ((b[3] <<  0) & 0x000000ff);
}

unsigned int RawData::Decode32bitWord(unsigned int word32bit = 0)
{
  //check data format
  unsigned int frame = word32bit & 0x80808080;
  if(frame != 0x80000000){
    std::cerr << __FILE__ << " " << __FUNCTION__ << " Frame Error! " << std::endl;
    std::cerr << "32 bit word: " << std::hex << word32bit << std::dec << std::endl;
    return 0;
  }

  return ((word32bit & 0x7f000000) >> 3) | 
         ((word32bit & 0x007f0000) >> 2) |
         ((word32bit & 0x00007f00) >> 1) |
         ((word32bit & 0x0000007f) >> 0);
}


//ADC High Gain
bool RawData::isAdcHg(unsigned int data = 0 )
{
    return (data & 0x00680000) == 0x00000000;
}


//ADC Low Gain
bool RawData::isAdcLg(unsigned int data = 0 )
{
    return (data & 0x00680000) == 0x00080000;
}

//TDC Leading
bool RawData::isTdcLeading(unsigned int data = 0 )
{
    return (data & 0x00601000) == 0x00201000;
}

//TDC Trailing
bool RawData::isTdcTrailing(unsigned int data = 0)
{
    return (data & 0x00601000) == 0x00200000;
}

//scaler function is not impletented yet in firmware (?) 
bool RawData::isScaler(unsigned int data = 0)
{
    return (data & 0x00600000) == 0x00400000;
}
