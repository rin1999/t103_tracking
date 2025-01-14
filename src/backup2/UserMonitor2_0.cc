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
  //DRS4
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
  int l1_tdc, l1_tdc1;

  std::vector<int> dclayer;
  int dcnhlayer;
  std::vector<int> dchitpat;
  int dcnhits[NumOfLayersDC];
  int dcwire[NumOfLayersDC];

  std::vector<std::vector<double> > dcltdc;
  std::vector<std::vector<double> > dcltdc_1st;
  std::vector<std::vector<double> > dcttdc;
  std::vector<std::vector<double> > dcttdc_1st;
  std::vector<std::vector<double> > dcwidth;
  std::vector<std::vector<double> > dcnhltdc;
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
  
  //**************************************************************************
  //DRS4
  {
#if check1
    std::cout << "DRS4***********************" << std::endl;
#endif

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
	  double adc = hit -> GetAmp(ia);
	  double amp = hit -> GetAmp(ia);
	  event.adc[segId-1].push_back(adc);
	  event.amp[segId-1].push_back(amp);
	}
	//TDC
	int sizeDt = hit->GetSize_Dt();
	for( int it=0; it<sizeDt; ++it ){
	  int dt = hit -> GetDt(it);
	  event.dt[segId-1].push_back(double(dt)*Tdc2Time);
	}
	int sizeWidth = hit->GetSize_Width();
	for( int it=0; it<sizeWidth; ++it ){
	  int width = hit -> GetWidth(it);
	  event.width[segId-1].push_back(double(width)*Tdc2Time);
	}
	event.l1_tdc = l1_tdc*Tdc2Time;
	event.l1_tdc1 = l1_tdc1*Tdc2Time;
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
    std::cout << "DC***********************" << std::endl;
#endif
    
    const int dctdclow = confMan->T0RangeLow();
    const int dctdchigh = confMan->T0RangeHigh();

    int nhlayerdc=0;
    int nhitsdc[NumOfLayersDC];
    
    bool tdcflagdc[NumOfLayersDC][MaxDCHits];
    
    for( int i=0; i<NumOfLayersDC; i++ ){
      nhitsdc[i]=0;
      for( int j=0; j<MaxWire; j++ ){
	tdcflagdc[i][j] = false;
      }
    }

    event.dcltdc.resize(MaxDCHits);
    event.dcltdc_1st.resize(MaxDCHits);
    event.dcttdc.resize(MaxDCHits);
    event.dcttdc_1st.resize(MaxDCHits);
    event.dcwidth.resize(MaxDCHits);
    event.dcnhltdc.resize(MaxDCHits);

    int nhits=0;
    int id=0;
    const HRTdcRHitContainer &cont =rawData->GetHRTRHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HRTdcRawHit *hit=cont[i];
      
      int detectorId = hit->DetectorId();
      int layerId = hit->LayerId()+1;
      int wireId = hit->SegmentId()+1;
      int UorD = hit->GetUorD();
	
      id=MaxWire*(layerId-1)+(wireId-1);

      int id2 = -1;
      if(layerId==1) id2=MaxWire*(layerId-1)+(wireId-1);
      if(layerId==2) id2=MaxWire*(layerId-1)+(wireId-1);
      if(layerId==3) id2=MaxWire*(layerId-1)+(wireId-1)-1;
      if(layerId==4) id2=id;
	
      if(id2<33){


	/*


	//TDC(leading)
	int sizelTdc = hit->GetSize_lTdc();	
	int tdc1st = -1;
	
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;
	  
	  if(id2%2==1) event.dcltdc[id].push_back(ltdc);
	  if(id2%2==0) event.dcttdc[id].push_back(ltdc);
 
	  if( dctdclow<ltdc && ltdc<dctdchigh ){
	    tdcflagdc[layerId-1][wireId-1] = true;
	  }
	  if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	  event.dcltdc_1st[id].push_back(tdc1st);
	}
	event.dcnhltdc[id].push_back(sizelTdc);
      
	if( tdcflagdc[layerId-1][wireId-1] ){
	  nhitsdc[layerId-1]++;
	  event.dclayer.push_back(layerId);
	  event.dcwire[layerId-1] = wireId;
	  
	  event.dchitpat.push_back(id2);
	}
	
	//TDC(trailing)
	int sizetTdc = hit->GetSize_tTdc();
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = double(hit -> GettTdc(itt))*Tdc2Time2;
	  
	  if(id2%2==0) event.dcltdc[id].push_back(ttdc);
	  if(id2%2==1) event.dcttdc[id].push_back(ttdc);

	  if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc < tdc1st ) tdc1st=ttdc;
	  event.dcttdc_1st[id].push_back(tdc1st);
	}


	*/


	//TDC
	int sizelTdc = hit->GetSize_lTdc();
	int sizetTdc = hit->GetSize_tTdc();	


	if(id2%2==1){

	    int ltdc1st = -1;
	    int ttdc1st = -1;

	//TDC(leading)
	  for( int irt=0; irt<sizelTdc; ++irt ){
	    int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;	


	    event.dcltdc[id].push_back(ltdc);
	    //if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	    if(ltdc < dctdchigh && ltdc > ltdc1st ) ltdc1st=ltdc;
	    //if( ltdc > ltdc1st ) ltdc1st=ltdc;
	    event.dcltdc_1st[id].push_back(ltdc1st);

	    if( dctdclow<ltdc && ltdc<dctdchigh ){
	      tdcflagdc[layerId-1][wireId-1] = true;
	    }
	  }

	//TDC(trailing)
	  for( int itt=0; itt<sizetTdc; ++itt ){
	    int ttdc = double(hit->GettTdc(itt))*Tdc2Time2;

	    event.dcttdc[id].push_back(ttdc);
	    //if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc < tdc1st ) tdc1st=ttdc;
	    if(ttdc < dctdchigh && ttdc > ttdc1st ) ttdc1st=ttdc;
	    //if( ttdc > ttdc1st ) ttdc1st=ttdc;

	    event.dcttdc_1st[id].push_back(ttdc1st);
	  }
	}
	  
	  
	if(id2%2==0){
	  
	  int ltdc1st = -1;
	  int ttdc1st = -1;

	//TDC(leading)
	  for( int itt=0; itt<sizetTdc; ++itt ){
	    int ttdc = double(hit->GettTdc(itt))*Tdc2Time2;


	    event.dcltdc[id].push_back(ttdc);
	    //if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	    if(ttdc < dctdchigh && ttdc > ltdc1st ) ltdc1st=ttdc;
	    //if( ttdc > ltdc1st ) ltdc1st=ttdc;
	    event.dcltdc_1st[id].push_back(ltdc1st);
	  }

	//TDC(trailing)
	  for( int irt=0; irt<sizelTdc; ++irt ){
	    int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;
	

	    event.dcttdc[id].push_back(ltdc);	  
	    //if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc < tdc1st ) tdc1st=ltdc;
	    if(ltdc < dctdchigh && ltdc > ttdc1st ) ttdc1st=ltdc;
	    //if( ltdc > ttdc1st ) ttdc1st=ltdc;
	    event.dcttdc_1st[id].push_back(ttdc1st);
 
	    if( dctdclow<ltdc && ltdc<dctdchigh ){
	      tdcflagdc[layerId-1][wireId-1] = true;
	    }
	  }
	}

	event.dcnhltdc[id].push_back(sizelTdc);
	if( tdcflagdc[layerId-1][wireId-1] ){
	  nhitsdc[layerId-1]++;
	  event.dclayer.push_back(layerId);
	  event.dcwire[layerId-1] = wireId;
	  
	  event.dchitpat.push_back(id2);
	}

	
	//ToT(width)
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;
	  for( int itt=0; itt<sizetTdc; ++itt ){
	    int ttdc = double(hit -> GettTdc(itt))*Tdc2Time2;
	    int twidth = ltdc-ttdc;
	    event.dcwidth[id].push_back(twidth);
	  }
	}
      }

      else {

	//TDC(leading)
	int sizelTdc = hit->GetSize_lTdc();	
	int tdc1st = -1;
	
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;
	  
	  event.dcltdc[id].push_back(ltdc);
 
	  if( dctdclow<ltdc && ltdc<dctdchigh ){
	    tdcflagdc[layerId-1][wireId-1] = true;
	  }
	  if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	  event.dcltdc_1st[id].push_back(tdc1st);
	}
	event.dcnhltdc[id].push_back(sizelTdc);
      
	if( tdcflagdc[layerId-1][wireId-1] ){
	  nhitsdc[layerId-1]++;
	  event.dclayer.push_back(layerId);
	  event.dcwire[layerId-1] = wireId;
	  
	  event.dchitpat.push_back(id2);
	}
	
	//TDC(trailing)
	int sizetTdc = hit->GetSize_tTdc();
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = double(hit -> GettTdc(itt))*Tdc2Time2;
	  
	  event.dcttdc[id].push_back(ttdc);
	}

	//ToT(width)
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = double(hit->GetlTdc(irt))*Tdc2Time2;
	  for( int itt=0; itt<sizetTdc; ++itt ){
	    int ttdc = double(hit -> GettTdc(itt))*Tdc2Time2;
	    int twidth = ltdc-ttdc;
	    event.dcwidth[id].push_back(twidth);
	  }
	}
      }

#if check1
	std::cout << "Layer = " << layerId 
		  << " Wire = " << wireId 
		  << std::endl;
#endif
    }

    for( int layer=0; layer<NumOfLayersDC; layer++ ){
      if( layer==0 ){
	if( tdcflagdc[0][0] || tdcflagdc[0][1] || tdcflagdc[0][2] ||
	    tdcflagdc[0][3] || tdcflagdc[0][4] || tdcflagdc[0][5] ||
	    tdcflagdc[0][6] || tdcflagdc[0][7] || tdcflagdc[0][8] ||
	    tdcflagdc[0][9] || tdcflagdc[0][10] ){
	  nhlayerdc++;
	}
      }
      if( layer==1 ){
	if( tdcflagdc[1][0] || tdcflagdc[1][1] || tdcflagdc[1][2] ||
	    tdcflagdc[1][3] || tdcflagdc[1][4] || tdcflagdc[1][5] ||
	    tdcflagdc[1][6] || tdcflagdc[1][7] || tdcflagdc[1][8] ||
	    tdcflagdc[1][9] ){
	  nhlayerdc++;
	}
      }
      if( layer==2 ){
	if( tdcflagdc[2][0] || tdcflagdc[2][1] || tdcflagdc[2][2] ||
	    tdcflagdc[2][3] || tdcflagdc[2][4] || tdcflagdc[2][5] ||
	    tdcflagdc[2][6] || tdcflagdc[2][7] || tdcflagdc[2][8] ||
	    tdcflagdc[2][9] || tdcflagdc[2][10] ){
	  nhlayerdc++;
	}
      }
    }
    event.dcnhlayer = nhlayerdc;

    for( int i=0; i<NumOfLayersDC; i++ ){
      event.dcnhits[i]=nhitsdc[i];
    }
  }

  tree->Fill();

  return true;
}

void EventMonitor2::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  //DRS4
  if( confMan->T0SimpleAna() ){
    for(int it=0; it<event.adc.size(); it++){
      event.adc[it].clear();
      event.amp[it].clear();
    }
    for(int it=0; it<event.tdc.size(); it++){
      event.dt[it].clear();
      event.width[it].clear();
    }
    event.l1_tdc  =-1;
    event.l1_tdc1 =-1;
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

  //DC
  event.dclayer.clear();
  event.dcnhlayer=-1;
  event.dchitpat.clear();
  
  for(int it=0; it<NumOfLayersDC; it++){
    event.dcnhits[it] =-1;
    event.dcwire[it] =-1;
  }
  
  event.dcltdc.clear();
  event.dcltdc_1st.clear();
  event.dcttdc.clear();
  event.dcttdc_1st.clear();
  event.dcwidth.clear();
  event.dcnhltdc.clear();
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

  //DRS4
  if( confMan->T0SimpleAna() ){
    tree->Branch("adc", &event.adc);
    tree->Branch("amp", &event.amp);
    tree->Branch("dt", &event.dt);
    tree->Branch("width", &event.width);
    tree->Branch("l1_tdc", &event.l1_tdc, "l1_tdc/I");
    tree->Branch("l1_tdc1", &event.l1_tdc1, "l1_tdc1/I");
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

  //DC
  tree->Branch("dclayer", &event.dclayer);
  tree->Branch("dcnhlayer", &event.dcnhlayer, "dcnhlayer/I");
  tree->Branch("dchitpat", &event.dchitpat);
  tree->Branch("dcnhits",  event.dcnhits, "dcnhits[8]/I");
  tree->Branch("dcwire",  event.dcwire, "dcwire[8]/I");
  
  tree->Branch("dcltdc", &event.dcltdc);
  tree->Branch("dcltdc_1st", &event.dcltdc_1st);
  tree->Branch("dcttdc", &event.dcttdc);
  tree->Branch("dcttdc_1st", &event.dcttdc_1st);
  tree->Branch("dcwidth", &event.dcwidth);
  tree->Branch("dcnhltdc", &event.dcnhltdc);

  return true;
}
