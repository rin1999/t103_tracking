/*
  UserMonitor2.cc

  2019/08 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "Decoder.hh"
#include "RootHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoRawHit.hh"

#include "TFile.h"
#include "TTree.h"

#define check1 0

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventMonitor2 : public VEvent
{
public:
  EventMonitor2();
  ~EventMonitor2();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( Decoder& gDec, int evnum );
  void InitializeEvent();

private:
  RawData *rawData;

};

EventMonitor2::EventMonitor2()
  : VEvent(),
    rawData(0)
{
}

EventMonitor2::~EventMonitor2()
{
  if (rawData) delete rawData;
}

struct Event{
  //HR-TDC
  int nh_dcltdc[NumOfHRTDC], nh_dcttdc[NumOfHRTDC];

  std::vector<std::vector<double> > dcltdc;
  std::vector<std::vector<double> > dcttdc;
  std::vector<std::vector<double> > width;

  int ev;

};
static Event event;

bool EventMonitor2::ProcessingBegin()
{
 return true;
}

bool EventMonitor2::ProcessingNormal( Decoder& gDec, int evnum )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits( gDec ) ) return false;
  //std::cout << "***" << std::endl;

  ConfMan *confMan = ConfMan::GetConfManager();

  //**************************************************************************
  //******************RawData
  TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();
  
  //**************************************************************************
  //HR-TDC
  {
#if check1
    std::cout << "HR-TDC***********************" << std::endl;
#endif
    event.dcltdc.resize(NumOfHRTDC);
    event.dcttdc.resize(NumOfHRTDC);
    event.width.resize(NumOfHRTDC);

    const HodoRHitContainer &cont =rawData->GetT0RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];

      int layerId = hit->LayerId()+1;
      int segId = hit->SegmentId()+1;

      //TDC leading
      int sizelTdc = hit->GetSize_lTdc();
      event.nh_dcltdc[segId-1] = sizelTdc;
      for( int it=0; it<sizelTdc; ++it ){
	int ltdc = (hit -> GetlTdc(it));

	if(segId%2==1) event.dcltdc[segId-1].push_back(double(ltdc)*Tdc2Time);
	if(segId%2==0) event.dcttdc[segId-1].push_back(double(ltdc)*Tdc2Time);
      }
      //TDC trailing
      int sizetTdc = hit->GetSize_tTdc();
      event.nh_dcttdc[segId-1] = sizetTdc;
      for( int it=0; it<sizetTdc; ++it ){
	int ttdc = (hit -> GettTdc(it));

	if(segId%2==0) event.dcltdc[segId-1].push_back(double(ttdc)*Tdc2Time);
	if(segId%2==1) event.dcttdc[segId-1].push_back(double(ttdc)*Tdc2Time);
      }

      //ToT(width)
      // int sizeWidth = hit->GetSize_Width();
      // for( int it=0; it<sizeWidth; ++it ){
      // 	int twidth = (hit -> GetWidth(it));
      // 	event.width[segId-1].push_back(double(twidth)*Tdc2Time);
      // }

#if check1
      std::cout << "Layer = " << layerId 
      		<< " Seg = " << segId 
      		<< std::endl;
#endif

    }
  }


  event.ev = evnum;
  tree->Fill();

  return true;
}

void EventMonitor2::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  for(int it=0; it<NumOfHRTDC; it++){
    event.nh_dcltdc[it] = -1;
    event.nh_dcttdc[it] = -1;
  }

  for(int it=0; it<event.dcltdc.size(); it++){
    event.dcltdc[it].clear();
  }
  for(int it=0; it<event.dcttdc.size(); it++){
    event.dcttdc[it].clear();
  }

  event.ev = 0;
}

bool EventMonitor2::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor2;
}

bool ConfMan:: InitializeHistograms()
{  
  ConfMan *confMan = ConfMan::GetConfManager();

  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  tree->Branch("nh_dcltdc",  event.nh_dcltdc, "nh_dcltdc[64]/I");
  tree->Branch("nh_dcttdc",  event.nh_dcttdc, "nh_dcttdc[64]/I");
  
  tree->Branch("dcltdc", &event.dcltdc);
  tree->Branch("dcttdc", &event.dcttdc);
  tree->Branch("width", &event.width);

  tree->Branch("event", &event.ev);
  
  return true;
}
