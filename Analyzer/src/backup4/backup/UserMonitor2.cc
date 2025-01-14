/*
  UserMonitor2.cc

  2018/10 K.Shirotori 
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
#include "HRTdcRawHit.hh"

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
  bool ProcessingNormal( Decoder& gDec );
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
  //T0
  int evDRS;

  std::vector<std::vector<double> > WaveForm;

  std::vector<std::vector<double> > adc;
  std::vector<std::vector<double> > amp;
  std::vector<std::vector<double> > bl;
  std::vector<std::vector<double> > peakx;

  std::vector<std::vector<double> > tdc;
  std::vector<std::vector<double> > dt;
  std::vector<std::vector<double> > tdc_2nd;
  std::vector<std::vector<double> > dt_2nd;
  std::vector<std::vector<double> > width;

  //HR-TDC
  int evHRT;

  std::vector<std::vector<double> > ltdc;
  std::vector<std::vector<double> > ttdc;

  int l1_tdc, l1_tdc1;

  double adc1[MaxHits], adc2[MaxHits], adc3[MaxHits], adc4[MaxHits];
  double amp1[MaxHits], amp2[MaxHits], amp3[MaxHits], amp4[MaxHits];
  double tdc1[MaxHits], tdc2[MaxHits], tdc3[MaxHits], tdc4[MaxHits];
  double tdct1[MaxHits], tdct2[MaxHits], tdct3[MaxHits], tdct4[MaxHits];
  double width1[MaxHits], width2[MaxHits], width3[MaxHits], width4[MaxHits];
  double ltdc1[MaxHits], ltdc2[MaxHits], ltdc3[MaxHits], ltdc4[MaxHits];
  double ttdc1[MaxHits], ttdc2[MaxHits], ttdc3[MaxHits], ttdc4[MaxHits];
  double lrftdc[MaxHits],trftdc[MaxHits];
};
static Event event;

bool EventMonitor2::ProcessingBegin()
{
 return true;
}

bool EventMonitor2::ProcessingNormal( Decoder& gDec )
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
  
  double time1=0., time2=0., diff=0.;

  int chx1=1, chx2=2, chs1=3, chs2=4;
  int chx1t=1, chx2t=3, chs1t=5, chs2t=7, chrft=16;

  //**************************************************************************
  //DRS4
  {
#if check1
    std::cout << "DRS4 ***********************" << std::endl;
#endif

    const int t0tdclow = confMan->T0RangeLow();
    const int t0tdchigh = confMan->T0RangeHigh();

    event.WaveForm.resize(NumOfCell);

    event.adc.resize(NumOfSegT0);
    event.amp.resize(NumOfSegT0);
    event.bl.resize(NumOfSegT0);
    event.peakx.resize(NumOfSegT0);
    
    event.tdc.resize(NumOfSegT0);
    event.dt.resize(NumOfSegT0);
    event.tdc_2nd.resize(NumOfSegT0);
    event.dt_2nd.resize(NumOfSegT0);
    event.width.resize(NumOfSegT0);

    int a1=0, a2=0, a3=0, a4=0;
    int t1=0, t2=0, t3=0, t4=0;
    int tt1=0, tt2=0, tt3=0, tt4=0;
    int w1=0, w2=0, w3=0, w4=0;

    const HodoRHitContainer &cont =rawData->GetT0RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];

      int detectorId = hit->DetectorId();
      int layerId = hit->LayerId()+1;
      int segId = hit->SegmentId()+1;
      int UorD = hit->GetUorD();
      int l1_tdc = hit -> GetL1_Tdc();
      int l1_tdc1 = hit -> GetL1_Tdc1();

      int eventNum = hit->GetEvNum();
      // std::cout<< "DRS4 Ev Num = " << eventNum << std::endl;
      event.evDRS = eventNum;

      if( confMan->T0SimpleAna() ){
	//ADC
	int sizeAdc = hit->GetSize_Amp();
	for( int ia=0; ia<sizeAdc; ++ia ){
      	  double adc = hit -> GetAdc(ia);
	  double amp = hit -> GetAmp(ia);
	  event.amp[segId-1].push_back(amp);
	  if(segId==chx1){
	    event.adc1[a1]=adc;
	    event.amp1[a1]=amp;
	    a1++;
	  }
	  if(segId==chx2){
	    event.adc2[a2]=adc;
	    event.amp2[a2]=amp;
	    a2++;
	  }
	  if(segId==chs1){
	    event.adc3[a3]=adc;
	    event.amp3[a3]=amp;
	    a3++;
	  }
	  if(segId==chs2){
	    event.adc4[a4]=adc;
	    event.amp4[a4]=amp;
	    a4++;
	  }
	}
	//TDC
	int sizeDt = hit->GetSize_Dt();
	for( int it=0; it<sizeDt; ++it ){
	  int dt = hit -> GetDt(it);
	  event.dt[segId-1].push_back(double(dt)*Tdc2Time);
	  if(segId==chx1){
	    event.tdc1[t1]=double(dt)*Tdc2Time;
	    if(t1==0){
	      time1=double(dt)*Tdc2Time;
	    }
	    t1++;
	  }
	  if(segId==chx2){
	    event.tdc2[t2]=double(dt)*Tdc2Time;
	    t2++;
	  }
	  if(segId==chs1){
	    event.tdc3[t3]=double(dt)*Tdc2Time;
	    t3++;
	  }
	  if(segId==chs2){
	    event.tdc4[t4]=double(dt)*Tdc2Time;
	    t4++;
	  }
	}
      	int sizeDt_2nd = hit->GetSize_Dt_2nd();
      	for( int it=0; it<sizeDt_2nd; ++it ){
      	  int dt_2nd = hit -> GetDt_2nd(it);
	  if(segId==chx1){
	    event.tdct1[tt1]=double(dt_2nd)*Tdc2Time;
	    tt1++;
	  }
	  if(segId==chx2){
	    event.tdct2[tt2]=double(dt_2nd)*Tdc2Time;
	    tt2++;
	  }
	  if(segId==chs1){
	    event.tdct3[tt3]=double(dt_2nd)*Tdc2Time;
	    tt3++;
	  }
	  if(segId==chs2){
	    event.tdct4[tt4]=double(dt_2nd)*Tdc2Time;
	    tt4++;
	  }
      	}
	int sizeWidth = hit->GetSize_Width();	
	for( int it=0; it<sizeWidth; ++it ){
	  int width = hit -> GetWidth(it);
	  event.width[segId-1].push_back(double(width)*Tdc2Time);
	  if(segId==chx1){
	    event.width1[w1]=double(width)*Tdc2Time;
	    w1++;
	  }
	  if(segId==chx2){
	    event.width2[w2]=double(width)*Tdc2Time;
	    w2++;
	  }
	  if(segId==chs1){
	    event.width3[w3]=double(width)*Tdc2Time;
	    w3++;
	  }
	  if(segId==chs2){
	    event.width4[w4]=double(width)*Tdc2Time;
	    w4++;
	  }
	}
	event.l1_tdc = l1_tdc*Tdc2Time;
	event.l1_tdc1 = l1_tdc1*Tdc2Time;
      }
      else{
      	// //Waveform
      	// int sizeWf = hit->GetSize_Wf();
      	// for( int iwf=0; iwf<sizeWf; ++iwf ){
      	//   double wf = hit -> GetWaveform(iwf);
      	//   event.WaveForm[segId-1].push_back(wf);
      	// }

      	//ADC
      	int sizeAdc = hit->GetSize_Amp();
      	for( int ia=0; ia<sizeAdc; ++ia ){
      	  double adc = hit -> GetAdc(ia);
      	  double amp = hit -> GetAmp(ia);
      	  double bl = hit -> GetBl(ia);
      	  double peakx = hit -> GetPeakX(ia);
	  
      	  event.adc[segId-1].push_back(adc);
      	  event.amp[segId-1].push_back(amp);
      	  event.bl[segId-1].push_back(bl);
      	  event.peakx[segId-1].push_back(peakx);
      	}
	
      	//TDC
      	int sizeTdc = hit->GetSize_Tdc();
      	for( int it=0; it<sizeTdc; ++it ){
      	  int tdc = (hit -> GetTdc(it));
      	  event.tdc[segId-1].push_back(double(tdc)*Tdc2Time);
      	}
      	int sizeDt = hit->GetSize_Dt();
      	for( int it=0; it<sizeDt; ++it ){
      	  int dt = hit -> GetDt(it);
      	  event.dt[segId-1].push_back(double(dt)*Tdc2Time);
      	}
      	int sizeTdc_2nd = hit->GetSize_Tdc_2nd();
      	for( int it=0; it<sizeTdc_2nd; ++it ){
      	  int tdc_2nd = hit -> GetTdc_2nd(it);
      	  event.tdc_2nd[segId-1].push_back(double(tdc_2nd)*Tdc2Time);
      	}
      	int sizeDt_2nd = hit->GetSize_Dt_2nd();
      	for( int it=0; it<sizeDt_2nd; ++it ){
      	  int dt_2nd = hit -> GetDt_2nd(it);
      	  event.dt_2nd[segId-1].push_back(double(dt_2nd)*Tdc2Time);
      	}
      	int sizeWidth = hit->GetSize_Width();
      	for( int it=0; it<sizeWidth; ++it ){
      	  int width = hit -> GetWidth(it);
      	  event.width[segId-1].push_back(double(width)*Tdc2Time);
      	}
      	event.l1_tdc = l1_tdc*Tdc2Time;
      	event.l1_tdc1 = l1_tdc1*Tdc2Time;
      }
      
#if check1
      std::cout << "***DRS4***" 
		<< "Layer = " << layerId 
      		<< " Seg = " << segId 
      		<< std::endl;
#endif
    }
    // std::cout<< a1 << " " << a2<< " "  << a3 << " "  << a4 << std::endl;
    // std::cout<< t1 << " " << t2<< " "  << t3 << " "  << t4 << std::endl;
  }

  //**************************************************************************
  //HRTdc
  {
#if check1
    std::cout << "HR-TDC ***********************" << std::endl;
#endif
    
    event.ltdc.resize(NumOfHRTDC);
    event.ttdc.resize(NumOfHRTDC);

    int lt1=0, lt2=0, lt3=0, lt4=0;
    int tt1=0, tt2=0, tt3=0, tt4=0;
    int lrf=0, trf=0;
    const HRTdcRHitContainer &cont =rawData->GetHRTRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HRTdcRawHit *hit=cont[i];
      
      int detectorId = hit->DetectorId();
      int layerId = hit->LayerId()+1;
      int segId = hit->SegmentId()+1;
      int UorD = hit->GetUorD();
      
      int eventNum = hit->GetEvNum();
      // std::cout<< "MZN Ev Num = " << eventNum << std::endl;
      event.evHRT = eventNum;
      
      if( confMan->T0SimpleAna() ){
	//HR-TDC
	int sizelTdc = hit->GetSize_lTdc();	
	//std::cout << sizelTdc << std::endl;
	for( int it=0; it<sizelTdc; ++it ){
	  unsigned int ltdc = hit -> GetlTdc(it);
	  event.ltdc[segId-1].push_back(double(ltdc)*Tdc2Time2);
	  if(segId==chx1t){
	    event.ltdc1[lt1]=double(ltdc)*Tdc2Time2;
	    if(lt1==0){
	      time2=double(ltdc)*Tdc2Time2;
	    }
	    lt1++;
	  }
	  if(segId==chx2t){
	    event.ltdc2[lt2]=double(ltdc)*Tdc2Time2;
	    lt2++;
	  }
	  if(segId==chs1t){
	    event.ltdc3[lt3]=double(ltdc)*Tdc2Time2;
	    lt3++;
	  }
	  if(segId==chs2t){
	    event.ltdc4[lt4]=double(ltdc)*Tdc2Time2;
	    lt4++;
	  }
	  if(segId==chrft){
	    event.lrftdc[lrf]=double(ltdc)*Tdc2Time2;
	    lrf++;
	  }
	}
	int sizetTdc = hit->GetSize_tTdc();
	for( int it=0; it<sizetTdc; ++it ){
	  unsigned int ttdc = hit -> GettTdc(it);
	  event.ttdc[segId-1].push_back(double(ttdc)*Tdc2Time2);
	  if(segId==chx1t){
	    event.ttdc1[tt1]=double(ttdc)*Tdc2Time2;
	    tt1++;
	  }
	  if(segId==chx2t){
	    event.ttdc2[tt2]=double(ttdc)*Tdc2Time2;
	    tt2++;
	  }
	  if(segId==chs1t){
	    event.ttdc3[tt3]=double(ttdc)*Tdc2Time2;
	    tt3++;
	  }
	  if(segId==chs2t){
	    event.ttdc4[tt4]=double(ttdc)*Tdc2Time2;
	    tt4++;
	  }
	  if(segId==chrft){
	    event.trftdc[trf]=double(ttdc)*Tdc2Time2;
	    trf++;
	  }
	}
      }
      else{
      	//HR-TDC
      	int sizelTdc = hit->GetSize_lTdc();
      	for( int it=0; it<sizelTdc; ++it ){
      	  unsigned int ltdc = hit -> GetlTdc(it);
      	  event.ltdc[segId-1].push_back(double(ltdc)*Tdc2Time2);
      	}
      	int sizetTdc = hit->GetSize_tTdc();
      	for( int it=0; it<sizetTdc; ++it ){
      	  unsigned int ttdc = hit -> GettTdc(it);
      	  event.ttdc[segId-1].push_back(double(ttdc)*Tdc2Time2);
      	}
      }
      
#if check1
      std::cout << "***HR-TDC***" 
		<< "Layer = " << layerId 
		<< " Seg = " << segId 
		<< std::endl;
#endif
    }
    // std::cout<< lt1 << " " << lt2<< " "  << lt3 << " "  << lt4 << std::endl;
    // std::cout<< tt1 << " " << tt2<< " "  << tt3 << " "  << tt4 << std::endl;
  }
  // diff=time1-time2;
  // if(time1>0 && time2>0){
  //   std::cout<< "****" << std::endl;
  //   std::cout<< "time1:time2 = " << time1 << ":"  << time2 << std::endl;
  //   std::cout<< "DT = " << diff << std::endl;
  // }

  // std::cout << "amp= " << event.amp.size() << std::endl;
  // std::cout << "tdc= " << event.ltdc.size() << std::endl;

  tree->Fill();

  return true;
}

void EventMonitor2::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  event.evDRS  =-1;
  event.evHRT  =-1;

  //DRS4
  if( confMan->T0SimpleAna() ){
    for(int it=0; it<event.adc.size(); it++){
      event.amp[it].clear();
    }
    for(int it=0; it<event.tdc.size(); it++){
      event.dt[it].clear();
      event.width[it].clear();
    }
    event.l1_tdc  =-1;
    event.l1_tdc1 =-1;

    for(int it=0; it<MaxHits; it++){
      event.adc1[it] = -999.0;
      event.adc2[it] = -999.0;
      event.adc3[it] = -999.0;
      event.adc4[it] = -999.0;
      event.amp1[it] = -999.0;
      event.amp2[it] = -999.0;
      event.amp3[it] = -999.0;
      event.amp4[it] = -999.0;
      event.tdc1[it] = -999.0;
      event.tdc2[it] = -999.0;
      event.tdc3[it] = -999.0;
      event.tdc4[it] = -999.0;
      event.tdct1[it] = -999.0;
      event.tdct2[it] = -999.0;
      event.tdct3[it] = -999.0;
      event.tdct4[it] = -999.0;
      event.width1[it] = -999.0;
      event.width2[it] = -999.0;
      event.width3[it] = -999.0;
      event.width4[it] = -999.0;
    }
  }
  else{
    // for(int it=0; it<event.WaveForm.size(); it++){
    //   event.WaveForm[it].clear();
    // }
    for(int it=0; it<event.adc.size(); it++){
      event.adc[it].clear();
      event.amp[it].clear();
      event.bl[it].clear();
      event.peakx[it].clear();
    }
    for(int it=0; it<event.tdc.size(); it++){
      event.tdc[it].clear();
      event.dt[it].clear();
      event.tdc_2nd[it].clear();
      event.dt_2nd[it].clear();
      event.width[it].clear();
    }
    event.l1_tdc  =-1;
    event.l1_tdc1 =-1;
  }

 //HRTdc
  if( confMan->T0SimpleAna() ){
    for(int it=0; it<event.ltdc.size(); it++){
      event.ltdc[it].clear();
    }
    for(int it=0; it<event.ttdc.size(); it++){
      event.ttdc.clear();
    }
    for(int it=0; it<MaxHits; it++){
      event.ltdc1[it] = -999.0;
      event.ltdc2[it] = -999.0;
      event.ltdc3[it] = -999.0;
      event.ltdc4[it] = -999.0;
      event.ttdc1[it] = -999.0;
      event.ttdc2[it] = -999.0;
      event.ttdc3[it] = -999.0;
      event.ttdc4[it] = -999.0;

      event.lrftdc[it] = -999.0;
      event.trftdc[it] = -999.0;
    }
  }
  else{
    for(int it=0; it<event.ltdc.size(); it++){
      event.ltdc[it].clear();
    }
    for(int it=0; it<event.ttdc.size(); it++){
      event.ttdc.clear();
    }
  }
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

  tree->Branch("evDRS", &event.evDRS, "evDRS/I");
  tree->Branch("evHRT", &event.evHRT, "evHRT/I");

  //DRS4
  if( confMan->T0SimpleAna() ){
    tree->Branch("amp", &event.amp);
    tree->Branch("dt", &event.dt);
    tree->Branch("width", &event.width);
    tree->Branch("l1_tdc", &event.l1_tdc, "l1_tdc/I");
    tree->Branch("l1_tdc1", &event.l1_tdc1, "l1_tdc1/I");

    tree->Branch("adc1", event.adc1, "adc1[16]/D");
    tree->Branch("adc2", event.adc2, "adc2[16]/D");
    tree->Branch("adc3", event.adc3, "adc3[16]/D");
    tree->Branch("adc4", event.adc4, "adc4[16]/D");
    tree->Branch("amp1", event.amp1, "amp1[16]/D");
    tree->Branch("amp2", event.amp2, "amp2[16]/D");
    tree->Branch("amp3", event.amp3, "amp3[16]/D");
    tree->Branch("amp4", event.amp4, "amp4[16]/D");
    tree->Branch("tdc1", event.tdc1, "tdc1[16]/D");
    tree->Branch("tdc2", event.tdc2, "tdc2[16]/D");
    tree->Branch("tdc3", event.tdc3, "tdc3[16]/D");
    tree->Branch("tdc4", event.tdc4, "tdc4[16]/D");
    tree->Branch("tdct1", event.tdct1, "tdct1[16]/D");
    tree->Branch("tdct2", event.tdct2, "tdct2[16]/D");
    tree->Branch("tdct3", event.tdct3, "tdct3[16]/D");
    tree->Branch("tdct4", event.tdct4, "tdct4[16]/D");
    tree->Branch("width1", event.width1, "width1[16]/D");
    tree->Branch("width2", event.width2, "width2[16]/D");
    tree->Branch("width3", event.width3, "width3[16]/D");
    tree->Branch("width4", event.width4, "width4[16]/D");
  }
  else{
    tree->Branch("WaveForm", &event.WaveForm );

    tree->Branch("adc", &event.adc);
    tree->Branch("amp", &event.amp);
    tree->Branch("bl", &event.bl);
    tree->Branch("peakx", &event.peakx);
    
    tree->Branch("tdc", &event.tdc);
    tree->Branch("dt", &event.dt);
    tree->Branch("tdc_2nd", &event.tdc_2nd);
    tree->Branch("dt_2nd", &event.dt_2nd);
    tree->Branch("width", &event.width);
    
    tree->Branch("l1_tdc", &event.l1_tdc, "l1_tdc/I");
    tree->Branch("l1_tdc1", &event.l1_tdc1, "l1_tdc1/I");
  }

  //HRTdc
  if( confMan->T0SimpleAna() ){
    tree->Branch("ltdc", &event.ltdc);
    tree->Branch("ttdc", &event.ttdc);  

    tree->Branch("ltdc1", event.ltdc1, "ltdc1[16]/D");
    tree->Branch("ltdc2", event.ltdc2, "ltdc2[16]/D");
    tree->Branch("ltdc3", event.ltdc3, "ltdc3[16]/D");
    tree->Branch("ltdc4", event.ltdc4, "ltdc4[16]/D");
    tree->Branch("ttdc1", event.ttdc1, "ttdc1[16]/D");
    tree->Branch("ttdc2", event.ttdc2, "ttdc2[16]/D");
    tree->Branch("ttdc3", event.ttdc3, "ttdc3[16]/D");
    tree->Branch("ttdc4", event.ttdc4, "ttdc4[16]/D");

    tree->Branch("lrftdc", event.lrftdc, "lrftdc[16]/D");
    tree->Branch("trftdc", event.trftdc, "trftdc[16]/D");
  }
  else{
    tree->Branch("ltdc", &event.ltdc);
    tree->Branch("ttdc", &event.ttdc);
  }

  return true;
}
