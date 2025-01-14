/*
  UserMonitor.cc

  2021/11 K.Shirotori 
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
#include "TrRawHit.hh"
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

class EventMonitor : public VEvent
{
public:
  EventMonitor();
  ~EventMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( Decoder& gDec );
  void InitializeEvent();

private:
  RawData *rawData;

};

EventMonitor::EventMonitor()
  : VEvent(),
    rawData(0)
{
}

EventMonitor::~EventMonitor()
{
  if (rawData) delete rawData;
}

struct Event{
  //SFT
  std::vector<int> layer;
  int nhlayer;
  int nhits[NumOfLayersSFT];
  std::vector<std::vector<double> > seg;

  std::vector<std::vector<double> > adch;
  std::vector<std::vector<double> > adcl;
  std::vector<std::vector<double> > tdcl;
  std::vector<std::vector<double> > tdct;
  std::vector<std::vector<double> > width;
  std::vector<std::vector<double> > nhtdcl;

  //T0
  std::vector<std::vector<double> > WaveForm;

  std::vector<std::vector<double> > adc;
  std::vector<std::vector<double> > amp;
  std::vector<std::vector<double> > bl;
  std::vector<std::vector<double> > peakx;

  // std::vector<std::vector<double> > tdc;
  // std::vector<std::vector<double> > dt;
  // std::vector<std::vector<double> > tdc_2nd;
  // std::vector<std::vector<double> > dt_2nd;
  // std::vector<std::vector<double> > width;

  // int l1_tdc, l1_tdc1;

  //int nhits; 
  int Nhits; 
  std::vector<int> ch;
  //  int ch;
  std::vector<std::vector<double> > ltdc;
  std::vector<std::vector<double> > ttdc;

  double ltdc1[64];
  double ttdc1[64];

};
static Event event;

bool EventMonitor::ProcessingBegin()
{
 return true;
}

bool EventMonitor::ProcessingNormal( Decoder& gDec )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits( gDec ) ){
    std::cerr << " Decoding error! " << std::endl;
    std::cerr << "skipping this event.... " << std::endl;
    std::cerr << std::endl;
    
    return false;
  } 
  //std::cout << "***" << std::endl;

  ConfMan *confMan = ConfMan::GetConfManager();

  //**************************************************************************
  //******************RawData
  TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  
  //**************************************************************************
  //SFT
  const int tdclow = confMan->TRangeLow();
  const int tdchigh = confMan->TRangeHigh();

  int nhlayer=0;
  int nhits[NumOfLayersSFT];

  bool tdcflag[NumOfLayersSFT][MaxFiber];
  for( int i=0; i<NumOfLayersSFT; i++ ){
    nhits[i]=0;
    for( int j=0; j<MaxFiber; j++ ){
      tdcflag[i][j]= false;
    }
  }
  
  event.seg.resize(NumOfLayersSFT);
  event.adch.resize(MaxHits);
  event.adcl.resize(MaxHits);
  event.tdcl.resize(MaxHits);
  event.tdct.resize(MaxHits);
  event.width.resize(MaxHits);
  event.nhtdcl.resize(MaxHits);
  
#if check1
  std::cout << "SFT***********************" << std::endl;
#endif
  
  int id=0;
  for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
    const TrRHitContainer &cont =rawData->GetSFTRHC(layer);
    int nh=cont.size();
    
    for( int i=0; i<nh; ++i ){
      TrRawHit *hit=cont[i];
      int layerId = hit->LayerId()+1;
      int fiberId = hit->FiberId()+1;
      int UorD = hit->GetUorD();
      
      id=MaxFiber*(layer-1)+(fiberId-1);
      
      //ADC
      int adch = hit->GetAdcH();
      if( UorD==0 ) event.adch[id].push_back(adch);


      int adcl = hit->GetAdcL();
      if( UorD==0 ) event.adcl[id].push_back(adcl);
      
	//TDC(leading)
      int sizelTdc = hit->GetSize_lTdc();
      for( int irt=0; irt<sizelTdc; ++irt ){
	int ltdc = hit->GetlTdc(irt);
	
	if( UorD==0 ){
	  event.tdcl[id].push_back(ltdc);
	  if( tdclow<ltdc && ltdc<tdchigh ){
	    tdcflag[layer-1][fiberId-1] = true;
	  }
	}
      }

      if( UorD==0 ) event.nhtdcl[id].push_back(sizelTdc);
      
      //if( tdcflag[layer-1][fiberId-1] && !(fiberId-1==15) ){
      if( tdcflag[layer-1][fiberId-1] ){
	nhits[layer-1]++;
	event.layer.push_back(layer);
	event.seg[layer-1].push_back(fiberId-1);
      }
      
	//TDC(trailing)
      int sizetTdc = hit->GetSize_tTdc();
      for( int itt=0; itt<sizetTdc; ++itt ){
	int ttdc = hit -> GettTdc(itt);
	if( UorD==0 ){
	  event.tdct[id].push_back(ttdc);
	}
      }
      
      //ToT(width)
      for( int irt=0; irt<sizelTdc; ++irt ){
	int ltdc = hit->GetlTdc(irt);
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = hit -> GettTdc(itt);
	  int twidth = ltdc-ttdc;
	  
	  if( UorD==0 ){
	    event.width[id].push_back(twidth);
	  }
	}
      }
      
#if check1
      std::cout << "Layer = " << layerId 
		<< " Fiber = " << fiberId 
		<< " UorD= " << UorD 
		<< " ADC= " << adch
		<< std::endl;
#endif
    }
  }
  
  bool layerflag0 = false;
  bool layerflag1 = false;
  bool layerflag2 = false; 
  for( int layer=0; layer<NumOfLayersSFT; layer++ ){
    for( int  fiber=0; fiber<NumOfFiber; fiber++ ){
      if ( layer==0 ){
	if( tdcflag[layer][fiber] ){
	  layerflag0 = true;
	}
      }
      if ( layer==1 ){
	if( tdcflag[layer][fiber] ){
	  layerflag1 = true;
	}
      }
      if ( layer==2 ){
	if( tdcflag[layer][fiber] ){
	  layerflag2 = true;
	}
      }
    }
  }
  if( layerflag0 ) nhlayer++;
  if( layerflag1 ) nhlayer++;
  if( layerflag2 ) nhlayer++;

  event.nhlayer = nhlayer;
  for( int i=0; i<NumOfLayersSFT; i++ ){
    event.nhits[i]=nhits[i];
  }


  //**************************************************************************
  //DRS4
  {
#if check1
    std::cout << "DRS4***********************" << std::endl;
#endif
    // const int t0tdclow = confMan->T0RangeLow();
    // const int t0tdchigh = confMan->T0RangeHigh();

    event.WaveForm.resize(NumOfCell);

    event.adc.resize(NumOfSegT0);
    event.amp.resize(NumOfSegT0);
    event.bl.resize(NumOfSegT0);
    event.peakx.resize(NumOfSegT0);
    
    // event.tdc.resize(NumOfSegT0);
    // event.dt.resize(NumOfSegT0);
    // event.tdc_2nd.resize(NumOfSegT0);
    // event.dt_2nd.resize(NumOfSegT0);
    // event.width.resize(NumOfSegT0);

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

      if( confMan->T0SimpleAna() ){
	//ADC
	int sizeAdc = hit->GetSize_Amp();
	for( int ia=0; ia<sizeAdc; ++ia ){
	  double amp = hit -> GetAmp(ia);
	  event.amp[segId-1].push_back(amp);
	}
	// //TDC
	// int sizeDt = hit->GetSize_Dt();
	// for( int it=0; it<sizeDt; ++it ){
	//   int dt = hit -> GetDt(it);
	//   event.dt[segId-1].push_back(double(dt)*Tdc2Time);
	// }
	// int sizeWidth = hit->GetSize_Width();
	// for( int it=0; it<sizeWidth; ++it ){
	//   int width = hit -> GetWidth(it);
	//   event.width[segId-1].push_back(double(width)*Tdc2Time);
	// }
	// event.l1_tdc = l1_tdc*Tdc2Time;
	// event.l1_tdc1 = l1_tdc1*Tdc2Time;
      }
      else{
	//Waveform
	int sizeWf = hit->GetSize_Wf();
	for( int iwf=0; iwf<sizeWf; ++iwf ){
	  double wf = hit -> GetWaveform(iwf);
	  event.WaveForm[segId-1].push_back(wf);
	}
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
	
	// //TDC
	// int sizeTdc = hit->GetSize_Tdc();
	// for( int it=0; it<sizeTdc; ++it ){
	//   int tdc = (hit -> GetTdc(it));
	//   event.tdc[segId-1].push_back(double(tdc)*Tdc2Time);
	// }
	// int sizeDt = hit->GetSize_Dt();
	// for( int it=0; it<sizeDt; ++it ){
	//   int dt = hit -> GetDt(it);
	//   event.dt[segId-1].push_back(double(dt)*Tdc2Time);
	// }
	// int sizeTdc_2nd = hit->GetSize_Tdc_2nd();
	// for( int it=0; it<sizeTdc_2nd; ++it ){
	//   int tdc_2nd = hit -> GetTdc_2nd(it);
	//   event.tdc_2nd[segId-1].push_back(double(tdc_2nd)*Tdc2Time);
	// }
	// int sizeDt_2nd = hit->GetSize_Dt_2nd();
	// for( int it=0; it<sizeDt_2nd; ++it ){
	//   int dt_2nd = hit -> GetDt_2nd(it);
	//   event.dt_2nd[segId-1].push_back(double(dt_2nd)*Tdc2Time);
	// }
	// int sizeWidth = hit->GetSize_Width();
	// for( int it=0; it<sizeWidth; ++it ){
	//   int width = hit -> GetWidth(it);
	//   event.width[segId-1].push_back(double(width)*Tdc2Time);
	// }
	// event.l1_tdc = l1_tdc*Tdc2Time;
	// event.l1_tdc1 = l1_tdc1*Tdc2Time;
      }
      
#if check1
      std::cout << "Layer = " << layerId 
		<< " Seg = " << segId 
		<< std::endl;
#endif
      
    }
  }

  //**************************************************************************
  //HRTdc
  {
#if check1
    std::cout << "HRTdc***********************" << std::endl;
#endif
    
    event.ch.resize(NumOfHRTDC);
    event.ltdc.resize(NumOfHRTDC);
    event.ttdc.resize(NumOfHRTDC);

    const int t0tdclow = confMan->T0RangeLow();
    const int t0tdchigh = confMan->T0RangeHigh();

    int Nhits=0;
    
    const HRTdcRHitContainer &cont =rawData->GetHRTRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HRTdcRawHit *hit=cont[i];

      int detectorId = hit->DetectorId();
      int layerId = hit->LayerId()+1;
      int segId = hit->SegmentId()+1;
      int UorD = hit->GetUorD();
      bool ltdcflag=false,ttdcflag=false;
      int ltdc1st =-1, ttdc1st = -1;

      if( confMan->T0SimpleAna() ){
	//HR-TDC
	int sizelTdc = hit->GetSize_lTdc();
	for( int it=0; it<sizelTdc; ++it ){
	  unsigned int ltdc = hit -> GetlTdc(it);
	  event.ltdc[segId-1].push_back(double(ltdc)*Tdc2Time2);
	  if( t0tdclow<(ltdc*Tdc2Time2)&&(ltdc*Tdc2Time2)<t0tdchigh ){
	    ltdcflag=true;
	    if( (ltdc*Tdc2Time2) > ltdc1st ){
	      ltdc1st=ltdc;
	    }
	  }
	}
	if(ltdcflag) event.ltdc1[segId-1]=double(ltdc1st)*Tdc2Time2;
	// std::cout<< "**LTDC" << std::endl;
	// std::cout<< segId-1 << ":"<<(double(ltdc1st)*Tdc2Time2) << std::endl;
	// std::cout<< "*" << std::endl;

	int sizetTdc = hit->GetSize_tTdc();
	for( int it=0; it<sizetTdc; ++it ){
	  unsigned int ttdc = hit -> GettTdc(it);
	  event.ttdc[segId-1].push_back(double(ttdc)*Tdc2Time2);
	  if( t0tdclow<(ttdc*Tdc2Time2)&&(ttdc*Tdc2Time2)<t0tdchigh ){
	    ttdcflag=true;
	    if( (ttdc*Tdc2Time2) > ttdc1st ){
	      ttdc1st=ttdc;
	    }
	  }
	}
      	if(ttdcflag) event.ttdc1[segId-1]=double(ttdc1st)*Tdc2Time2;

	// std::cout<< "**TTDC" << std::endl;
	// std::cout<< segId-1 << ":"<<(double(ttdc1st)*Tdc2Time2) << std::endl;
      }
      else{
	//HR-TDC
	int sizelTdc = hit->GetSize_lTdc();
	for( int it=0; it<sizelTdc; ++it ){
	  unsigned int ltdc = hit -> GetlTdc(it);
	  event.ltdc[segId-1].push_back(double(ltdc)*Tdc2Time2);
	  if( t0tdclow<(ltdc*Tdc2Time2)&&(ltdc*Tdc2Time2)<t0tdchigh ) ltdcflag=true;
	}
	int sizetTdc = hit->GetSize_tTdc();
	for( int it=0; it<sizetTdc; ++it ){
	  unsigned int ttdc = hit -> GettTdc(it);
	  event.ttdc[segId-1].push_back(double(ttdc)*Tdc2Time2);
	}
      }
    
      //      if( ltdcflag && segId<10 ){
      //      if( ltdcflag ){
      if( ltdcflag && segId<129){
	Nhits++;
	event.ch.push_back(segId-1);
      }

#if check1
      std::cout << "Layer = " << layerId 
		<< " Seg = " << segId 
		<< std::endl;
#endif
    }
    // event.nhits=nhits;
       event.Nhits=Nhits;

  }

  tree->Fill();

  return true;
}

void EventMonitor::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  //SFT
  event.layer.clear();
  event.seg.clear();
  
  event.nhlayer=-1;
  
  for(int it=0; it<NumOfLayersSFT; it++){
    event.nhits[it] =-1;
  }
  
  event.adch.clear();
  event.adcl.clear();
  event.tdcl.clear();
  event.tdct.clear();
  event.width.clear();
  event.nhtdcl.clear();

  //DRS4
  if( confMan->T0SimpleAna() ){
    for(int it=0; it<event.adc.size(); it++){
      event.amp[it].clear();
    }
    // for(int it=0; it<event.tdc.size(); it++){
    //   event.dt[it].clear();
    //   event.width[it].clear();
    // }
    // event.l1_tdc  =-1;
    // event.l1_tdc1 =-1;
  }
  else{
    for(int it=0; it<event.WaveForm.size(); it++){
      event.WaveForm[it].clear();
    }
    for(int it=0; it<event.adc.size(); it++){
      event.adc[it].clear();
      event.amp[it].clear();
      event.bl[it].clear();
      event.peakx[it].clear();
    }
    // for(int it=0; it<event.tdc.size(); it++){
    //   event.tdc[it].clear();
    //   event.dt[it].clear();
    //   event.tdc_2nd[it].clear();
    //   event.dt_2nd[it].clear();
    //   event.width[it].clear();
    // }
    // event.l1_tdc  =-1;
    // event.l1_tdc1 =-1;
  }

  //HRTdc
  if( confMan->T0SimpleAna() ){
    for(int it=0; it<event.ltdc.size(); it++){
      event.ltdc[it].clear();
    }
    for(int it=0; it<event.ttdc.size(); it++){
      event.ttdc.clear();
    }
    event.ch.clear();
    // event.nhits =-1;

    for(int it=0; it<64; it++){
      event.ltdc1[it] = -999.0;
      event.ttdc1[it] = -999.0;
    } 
  }
  else{
    for(int it=0; it<event.ltdc.size(); it++){
      event.ltdc[it].clear();
    }
    for(int it=0; it<event.ttdc.size(); it++){
      event.ttdc.clear();
    }
    event.ch.clear();
    // event.nhits =-1;
  }
}

bool EventMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor;
}

bool ConfMan:: InitializeHistograms()
{  
  ConfMan *confMan = ConfMan::GetConfManager();

  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //SFT
  tree->Branch("layer", &event.layer);
  tree->Branch("seg", &event.seg);
  
  tree->Branch("nhlayer", &event.nhlayer, "nhlayer/I");
  tree->Branch("nhits",  event.nhits, "nhits[12]/I");
  
  tree->Branch("adch", &event.adch);
  tree->Branch("adcl", &event.adcl);
  tree->Branch("tdcl", &event.tdcl);
  tree->Branch("tdct", &event.tdct);
  tree->Branch("width", &event.width);
  tree->Branch("nhtdcl", &event.nhtdcl);

  //DRS4
  if( confMan->T0SimpleAna() ){
    tree->Branch("amp", &event.amp);
    // tree->Branch("dt", &event.dt);
    // tree->Branch("width", &event.width);
    // tree->Branch("l1_tdc", &event.l1_tdc, "l1_tdc/I");
    // tree->Branch("l1_tdc1", &event.l1_tdc1, "l1_tdc1/I");
  }
  else{
    tree->Branch("WaveForm", &event.WaveForm );

    tree->Branch("adc", &event.adc);
    tree->Branch("amp", &event.amp);
    tree->Branch("bl", &event.bl);
    tree->Branch("peakx", &event.peakx);
    
    // tree->Branch("tdc", &event.tdc);
    // tree->Branch("dt", &event.dt);
    // tree->Branch("tdc_2nd", &event.tdc_2nd);
    // tree->Branch("dt_2nd", &event.dt_2nd);
    // tree->Branch("width", &event.width);
    
    // tree->Branch("l1_tdc", &event.l1_tdc, "l1_tdc/I");
    // tree->Branch("l1_tdc1", &event.l1_tdc1, "l1_tdc1/I");
  }

  //HRTdc
  if( confMan->T0SimpleAna() ){
    tree->Branch("ltdc", &event.ltdc);
    tree->Branch("ttdc", &event.ttdc); 
    tree->Branch("ch", &event.ch);
    // tree->Branch("nhits", &event.nhits, "nhits/I");
    tree->Branch("Nhits", &event.Nhits, "Nhits/I");

    tree->Branch("ltdc1", event.ltdc1, "ltdc1[64]/D");
    tree->Branch("ttdc1", event.ttdc1, "ttdc1[64]/D");
  }
  else{
    tree->Branch("ltdc", &event.ltdc);
    tree->Branch("ttdc", &event.ttdc);
    tree->Branch("ch", &event.ch);
    // tree->Branch("nhits", &event.nhits, "nhits/I");
  }

  return true;
}
