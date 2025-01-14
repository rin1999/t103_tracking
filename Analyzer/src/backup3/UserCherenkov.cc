/*
  UserCherenkov.cc
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>

#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
//______________________________________________________________________________
VEvent::VEvent()
{
}
//______________________________________________________________________________
VEvent::~VEvent()
{
}
//______________________________________________________________________________
class EventCherenkov
  :public VEvent
{
private:
  RawData* rawData;

public:
  EventCherenkov();
  ~EventCherenkov();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
};
//______________________________________________________________________________
EventCherenkov::EventCherenkov()
  :VEvent(),
   rawData(0)
{
}
//______________________________________________________________________________
EventCherenkov::~EventCherenkov()
{
  if(rawData)
    delete rawData;
}
//______________________________________________________________________________
bool 
EventCherenkov::ProcessingBegin()
{
  ConfMan::InitializeHistograms();
  return true;
}
//______________________________________________________________________________
bool
EventCherenkov::ProcessingNormal()
{
  rawData = new RawData;
  rawData ->DecodeHits();

  //GC
  {
    const ChereRHitContainer &cont=rawData->GetGCRawHC();
    int nh=cont.size();
    //    HF1(100000,double(nh));
    int nh1=0,nh2=0;
    for(int i=0;i<nh;++i)
      {
	ChereRawHit *Hit=cont[i];
	int adc=hit->GetAdc();
	int tdc=hit->GetTdc();	

	HF1(100102,double(tdc));
	HF1(100101,double(adc));
      }

  }
  return true;
}
//______________________________________________________________________________
bool
EventCherenkov::ProcessingEnd()
{
  gFile->Write();
  gFile->Close();
  return true;
}
//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new EventCherenkov;
}
//______________________________________________________________________________
const int NbinGCAdc=500;
const int NbinGCTdc=4096;

const int MinGCAdc=0;
const int MaxGCAdc=500;

const int MinGCTdc=0;
const int MaxGCTdc=4096;

//______________________________________________________________________________
bool ConfMan::InitializeHistograms()
{
  //GC
  HB1(100101,"GC ADC",NbinGCAdc,MinGCAdc,MaxGCAdc);
  HB1(100102,"GC TDC",NbinGCTdc,MinGCTdc,MaxGCTdc);
  return true;
}
